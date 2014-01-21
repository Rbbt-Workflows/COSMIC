require 'rbbt'
require 'rbbt/workflow'
require 'rbbt/knowledge_base'
require 'progress-bar'

Workflow.require_workflow "Genomics"
require 'rbbt/entity/genomic_mutation'
require 'rbbt/entity/mutated_isoform'
require 'rbbt/entity/protein'

Workflow.require_workflow "Appris"
require 'appris'


module COSMIC
  extend Workflow

  class << self
    attr_accessor :knowledge_base, :organism
  end

  self.organism = "Hsa/jan2013"
  self.knowledge_base = KnowledgeBase.new Rbbt.var.knowledge_base.COSMIC, self.organism

  COSMIC.claim COSMIC.gene_damage_analysis, :proc do
    Workflow.require_workflow "MutEval"
    require 'db_nsfp'
    require 'rsruby'

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", 
                    :fields => ["Avg. damage score", "Bg. Avg. damage score", "T-test p-value"], 
                    :type => :list, :cast => :to_f, :namespace => COSMIC.organism)

    damage_fields = DbNSFP.database.fields.select{|f| f =~ /converted/ }
    db = COSMIC.knowledge_base.get_database(:gene_principal_isoform_mutations)
    db.with_monitor :desc => "Damage analysis using DbNSFP", :step => 10000 do
      db.through do |gene,mis|
        next if mis.empty?
        protein = mis.first.protein

        dbNSFP_tsv = DbNSFP.database.get_prefix(protein).slice(damage_fields)

        all_damage_scores = dbNSFP_tsv.collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact
        damage_scores = dbNSFP_tsv.select(:key => mis).collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact

        if damage_scores.length < 3
          damage_score_pvalue = 1
        else
          damage_score_pvalue = RSRuby.instance.t_test(damage_scores, all_damage_scores,"greater")["p.value"]
        end

        tsv[gene] = [Misc.mean(all_damage_scores), Misc.mean(damage_scores), damage_score_pvalue]
      end
    end

    tsv.to_s
  end
end

require 'rbbt/sources/COSMIC'

COSMIC.knowledge_base.register :sample_mutations, COSMIC.mutations, :source => "Sample name=>Sample", :fields => ["Genomic Mutation"], :merge => true, :unnamed => true, :namespace => COSMIC.organism

COSMIC.knowledge_base.register :mutation_genes do
  all_mutations = COSMIC.knowledge_base.get_database(:sample_mutations).values.compact.flatten.uniq
  GenomicMutation.setup(all_mutations, "All COSMIC mutation", "Hsa/jan2013", false)

  tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => COSMIC.organism, :unnamed => true)
  all_mutations.each do |mutation|
    tsv[mutation] = mutation.genes
  end
  tsv
end

COSMIC.knowledge_base.register :mutation_isoforms do
  all_mutations = COSMIC.knowledge_base.get_database(:sample_mutations).values.compact.flatten.uniq
  GenomicMutation.setup(all_mutations, "All COSMIC mutation", COSMIC.organism, false)

  tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Mutated Isoforms"], :type => :flat, :namespace => COSMIC.organism)
  all_mutations.each do |mutation|
    tsv[mutation] = mutation.mutated_isoforms
  end
  tsv
end

COSMIC.knowledge_base.register :gene_principal_isoform_mutations do
  all_mutations = COSMIC.knowledge_base.get_database(:sample_mutations).values.compact.flatten.uniq
  organism = COSMIC.organism
  GenomicMutation.setup(all_mutations, "All COSMIC mutation", organism, false)

  tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Mutated Isoform"], :type => :flat, :namespace => organism)


  all_mis = Annotated.purge(COSMIC.knowledge_base.get_index(:mutation_isoforms).keys).collect{|p| p.partition("~").last }.flatten.uniq
  protein_genes = Organism.transcripts(organism).tsv :persist => true, :key_field => "Ensembl Protein ID", :fields => ["Ensembl Gene ID"], :type => :single
  mi2genes = Misc.process_to_hash(all_mis) do |mis|
    mis.collect do |mi|
      protein, sep, change = mi.partition ":"
      next unless protein =~ /ENSP/
        next unless Appris::PRINCIPAL_ISOFORMS.include? protein
      protein_genes[protein]
    end
  end

  bar = Progress::Bar.new all_mutations.length, 0, 10000, "Mutations"
  all_mutations.each do |mutation|
    bar.tick
    matches = COSMIC.knowledge_base.children(:mutation_isoforms, mutation)
    next if matches.nil? or matches.empty?
    matches.target.each do |mi|
      gene = mi2genes[mi]
      next if gene.nil?
      tsv[gene] ||= []
      tsv[gene] << mi
    end
  end

  tsv
end

COSMIC.gene_damage_analysis.produce if __FILE__ == $0

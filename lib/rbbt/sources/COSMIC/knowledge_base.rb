require 'rbbt/knowledge_base'

module COSMIC
  class << self 
    attr_accessor :knowledge_base
  end
  self.knowledge_base = KnowledgeBase.new Rbbt.var.knowledge_base.COSMIC, self.organism
end

COSMIC.knowledge_base.register :sample_mutations, COSMIC.mutations, 
  :source => "Sample name=>Sample", :fields => ["Genomic Mutation"], :merge => true, :unnamed => true, :namespace => COSMIC.organism

COSMIC.knowledge_base.register :mutation_genes do
  Workflow.require_workflow "Sequence"

  all_mutations = COSMIC.knowledge_base.get_database(:sample_mutations).values.compact.flatten.uniq
  mutation2genes = Sequence.job(:genes_at_genomic_positions, "COSMIC", :positions => all_mutations, :organism => COSMIC.organism).run

  tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => COSMIC.organism, :unnamed => true)
  all_mutations.each do |mutation|
    tsv[mutation] = mutation2genes[mutation]
  end
  tsv
end

COSMIC.knowledge_base.register :mutation_isoforms do
  Workflow.require_workflow "Sequence"

  all_mutations = COSMIC.knowledge_base.get_database(:sample_mutations).values.compact.flatten.uniq
  mutation2mis = Sequence.job(:mutated_isoforms_for_genomic_mutations, "COSMIC", :mutations => all_mutations, :organism => COSMIC.organism).run

  tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Mutated Isoforms"], :type => :flat, :namespace => COSMIC.organism)
  all_mutations.each do |mutation|
    tsv[mutation] = mutation2mis[mutation]
  end
  tsv
end

COSMIC.knowledge_base.register :gene_principal_isoform_mutations do
  Workflow.require_workflow "Appris"
  require 'rbbt/sources/appris'

  all_mutations = COSMIC.knowledge_base.get_database(:sample_mutations).values.compact.flatten.uniq
  organism = COSMIC.organism

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

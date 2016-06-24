require 'rbbt/knowledge_base'

module COSMIC
  class << self 
    attr_accessor :knowledge_base
  end
  self.knowledge_base = KnowledgeBase.new Rbbt.var.COSMIC.knowledge_base, self.organism
end

COSMIC.knowledge_base.register :sample_mutations, COSMIC.mutations, 
  :source => "Sample name=~Sample", :target => ["Genomic Mutation"], :fields => ["Mutation zygosity"], :merge => true, :unnamed => true, :namespace => COSMIC.organism, :type => :double

COSMIC.knowledge_base.register :mutation_genes do
  Workflow.require_workflow "Sequence"

  all_mutations = COSMIC.knowledge_base.get_database(:sample_mutations).values.compact.flatten.uniq
  mutation2genes = Sequence.job(:genes, "COSMIC", :positions => all_mutations, :organism => COSMIC.organism).run

  tsv = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"], :type => :flat, :namespace => COSMIC.organism, :unnamed => true)
  all_mutations.each do |mutation|
    tsv[mutation] = mutation2genes[mutation]
  end

  tsv
end

COSMIC.knowledge_base.register :mutation_isoforms do
  Workflow.require_workflow "Sequence"

  sample_mutations = COSMIC.knowledge_base.get_database(:sample_mutations)
  all_mutations = sample_mutations.values.compact.flatten.uniq
  mutation2mis = Sequence.job(:mutated_isoforms_fast, "COSMIC", :mutations => all_mutations, :organism => COSMIC.organism, :watson => false)
  mutation2mis.run(true).join

  mutation2mis.path
end

COSMIC.knowledge_base.register :mutation_protein_changes do
  database = COSMIC.knowledge_base.get_database(:mutation_isoforms)

  TSV.traverse database, :into => :dumper,
    :key_field => "Genomic Mutation", :fields => ["Ensembl Protein ID", "Change"],
    :namespace => COSMIC.organism, :type => :double do |mutation, mis|

    values = mis.flatten.collect do |mi|
      protein, _sep, change = mi.partition ":"
      next unless change =~ /[A-Z*]\d+(?:[A-Z*]|FrameShift)/
      [protein, change]
    end.compact


    [mutation, Misc.zip_fields(values)]
  end
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

  all_mutations.each do |mutation|
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

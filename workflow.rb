require 'rbbt'
require 'rbbt/workflow'

module COSMIC
  extend Workflow

  class << self
    attr_accessor :organism
  end

  self.organism = "Hsa/jan2013"

  helper :organism do
    COSMIC.organism
  end

  input :mutated_isoforms, :array, "Mutated Isoforms"
  input :genes, :array, "Genes of interest"
  input :soft_match, :boolean, "Consider samples with different mutations but in the same residue", false
  desc "Find samples where particular protein mutations co-occurr with mutations in given genes"
  task :coocurrence_matrix => :tsv do |mutated_isoforms,genes,soft_match|

    Workflow.require_workflow "Proteomics"
    # Translate genes
    ensg_index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :persist => true
    ensg2name = Organism.identifiers(organism).index :fields => ["Ensembl Gene ID"], :target => "Associated Gene Name", :persist => true
    genes = ensg_index.chunked_values_at(genes).sort

    # Find samples
    mutation_info = Proteomics.job(:annotate_mi, nil, :mutated_isoforms => mutated_isoforms, :database => "COSMIC").run
    samples = mutation_info.column("Sample name").values.flatten.uniq

    # Determing which genes are mutated in each sample
    ensp2ensg = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :persist => true, :fields => ["Ensembl Protein ID"]
    sample_mutations = COSMIC.mutations.tsv :key_field => "Sample name", :fields => ["Genomic Mutation"], :type => :flat, :persist => true, :merge => true, :monitor => true

    tsv = TSV.setup({}, :key_field => "Sample", :fields => ["Ensembl Gene ID", "Mutated Isoform"], :type => :double, :organism => organism)
    TSV.traverse samples, :into => tsv, :type => :array do |sample|
      next unless sample_mutations.include? sample
      genomic_mutations = sample_mutations[sample]
      next if genomic_mutations.nil?
      mis = Sequence.job(:mutated_isoforms_fast, nil, :mutations => genomic_mutations, :watson => false, :organism => organism, :non_synonymous => true).run.values.flatten.compact.uniq
      next unless soft_match or (mis & mutated_isoforms).any?
      proteins = mis.collect{|mi| mi.partition(":").first }

      hit_genes = ensp2ensg.chunked_values_at proteins

      hit_mis = mis.select{|mi| p = mi.partition(":").first; g = ensp2ensg[p]; genes.include? g}
      [sample, [hit_genes & genes, hit_mis]]
    end

    # Build incidence matrix and add PMIDs
    index2 = COSMIC.mutations.tsv :key_field => "Sample name", :fields => "Pubmed_PMID", :type => :flat, :persist => true, :merge => true, :monitor => true
    matrix = TSV.setup({}, :key_field => "Sample name", :fields => genes, :type => :double, :organism => organism)
    tsv.each do |sample,values|
      hit_genes, hit_mis = values
      #matrix[sample] = genes.collect{|g| [hit_genes.include?(g) ? 'true' : 'false']}
      gene_mis = {}
      hit_mis.each do |mi|
         p = mi.partition(":").first
         g = ensp2ensg[p]
         gene_mis[g] ||= []
         gene_mis[g] << mi
      end
      matrix[sample] = genes.collect{|g| gene_mis[g] }
    end

    matrix.fields = ensg2name.chunked_values_at matrix.fields

    res = matrix.attach index2
    res.fields[-1] = "PMID"
    res = res.process "PMID" do |value, key, values|
      value.flatten.uniq
    end
    res
  end

  export_synchronous :coocurrence_matrix
end

require 'rbbt/sources/COSMIC'
require 'rbbt/knowledge_base/COSMIC'


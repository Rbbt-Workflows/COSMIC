require 'rbbt'
require 'rbbt/workflow'

module COSMIC
  extend Workflow

  class << self
    attr_accessor :organism
  end

  self.organism = "Hsa/jan2013"
end

require 'rbbt/sources/COSMIC'

#Workflow.require_workflow "Genomics"
#require 'rbbt/entity/genomic_mutation'
#require 'rbbt/entity/mutated_isoform'

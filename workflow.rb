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


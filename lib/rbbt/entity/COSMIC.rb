module Gene
  property :COSMIC_rsids => :single2array do
    COSMIC.rsid_index(organism, chromosome)[self.chr_range]
  end

  property :COSMIC_mutations => :single2array do
    GenomicMutation.setup(COSMIC.mutation_index(organism).values_at(*self.COSMIC_rsids).uniq, "COSMIC mutations over #{self.name || self}", organism, false)
  end
end if defined? Entity and defined? Gene and Entity === Gene

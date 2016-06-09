require 'rbbt'
require 'rbbt/tsv'
require 'rbbt/resource'
require 'rbbt/util/misc/omics'

module COSMIC
  extend Resource
  self.subdir = "share/databases/COSMIC"

  def self.organism
    Organism.default_code "Hsa"
  end

  COSMIC.claim COSMIC.CosmicMutantExport, :proc do |filename|
    url = "sftp://sftp-cancer.sanger.ac.uk/cosmic/grch37/cosmic/v77/CosmicMutantExport.tsv.gz"
    raise "Follow #{ url } and place the file uncompressed in #{filename}"
  end

  COSMIC.claim COSMIC.CosmicResistanceMutations, :proc do |filename|
    url = "sftp://sftp-cancer.sanger.ac.uk/cosmic/grch37/cosmic/v77/CosmicResistanceMutations.tsv.gz"
    raise "Follow #{ url } and place the file uncompressed in #{filename}"
  end

  COSMIC.claim COSMIC.sample_info, :proc do |file|
    site_and_histology_fieds = TSV.parse_header(COSMIC.mutations).fields.select{|f| f =~ /site|histology/i }
    tsv = COSMIC.mutations.tsv(:key_field => "Sample name", :fields => site_and_histology_fieds, :type => :list)
    tsv.to_s
  end

  COSMIC.claim COSMIC.mutations, :proc do |directory|
    url = COSMIC.CosmicMutantExport.produce.find
    #stream = CMD.cmd('awk \'BEGIN{FS="\t"} { if ($8 != "" && $8 != "Mutation ID") { sub($8, "COSM" $8 ":" $3)}; print}\'', :in => Open.open(url), :pipe => true)

    #all_fields = TSV.parse_header(url, :header_hash => "").all_fields
    #mutation_id_field, sample_id_field = ["Mutation ID", "ID_sample"].collect{|f| all_fields.index(f) + 1}
    #stream = CMD.cmd("awk 'BEGIN{FS=\"\\t\"} { if ($#{mutation_id_field} != \"\" && $#{mutation_id_field} != \"Mutation ID\") { sub($#{mutation_id_field}, \"COSM:\" $#{mutation_id_field} \":\" $#{sample_id_field})}; print}'", :in => Open.open(url), :pipe => true)

    stream = Open.open(url)
    parser = TSV::Parser.new stream, :header_hash => "", :key_field => "Mutation ID", :type => :list

    dumper = TSV::Dumper.new parser.options.merge({:fields => parser.fields + ["Genomic Mutation"]})
    dumper.init

    pos_i = parser.identify_field "Mutation genome position"
    cds_i = parser.identify_field "Mutation CDS"
    sample_i = parser.identify_field "ID_sample"
    dumper = TSV.traverse parser, :type => :list, :into => dumper, :bar => url do |mid, values|
      position = values[pos_i]
      cds = values[cds_i]
      sample = values[sample_i]

      mid = ["COSM", mid, sample] * ":"

      if position.nil? or position.empty?
        nil
      else
        position = position.split("-").first

        chr, pos = position.split(":")
        chr = "X" if chr == "23"
        chr = "Y" if chr == "24"
        chr = "MT" if chr == "25"
        position = [chr, pos ] * ":"
        next if chr.length >= 3

        if cds.nil?
          next
        else
          change = Misc.translate_dna_mutation_hgvs2rbbt(cds)
        end
      end
      genomic_mutation = [chr,pos,change] * ":"
      [mid,values+[genomic_mutation]]
    end
    dumper.stream
  end

  COSMIC.claim COSMIC.mi_drug_resistance, :proc do |filename|
    tsv = COSMIC.CosmicResistanceMutations.tsv :key_field => "Transcript", :fields => ["AA Mutation", "Drug Name","Sample ID", "Pubmed Id", "Zygosity"], :type => :double, :merge => true, :header_hash => ''
    res = TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Drug name","Sample","PMID", "Zygosity"], :type => :double)
    organism = "Hsa/feb2014"
    enst2ensp = Organism.transcripts(organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Protein ID"], :type => :single, :merge => true, :persist => true
    TSV.traverse tsv, :into => res, :bar => true do |transcript, values|
      protein = enst2ensp[transcript]
      new = []
      Misc.zip_fields(values).each do |aa_mutation,drug,sample,pmid,zygosity|
        mutation = Misc.translate_prot_mutation_hgvs2rbbt(aa_mutation)
        mi = [protein,mutation] * ":"
        new << [mi, [drug,sample,pmid,zygosity]]
      end
      new.extend MultipleResult
      new
    end
    res.to_s
  end

  COSMIC.claim COSMIC.mutations_hg18, :proc do |filename|
    require 'rbbt/sources/organism'
    file = COSMIC.mutations.open
    begin

      while (line = file.gets) !~ /Genomic Mutation/; end
      fields = line[1..-2].split("\t")
      mutation_pos = fields.index "Genomic Mutation"

      mutations = CMD.cmd("grep -v '^#'|cut -f #{mutation_pos + 1}|sort -u", :in => COSMIC.mutations.open).read.split("\n").select{|m| m.include? ":" }

      translations = Misc.process_to_hash(mutations){|mutations| Organism.liftOver(mutations, self.organism, "Hsa/may2009")}

      File.open(filename, 'w') do |f|
        f.puts "#: :type=:list#:namespace=Hsa/may2009"
        f.puts "#" + fields * "\t"
        while line = file.gets do
          next if line[0] == "#"[0]
          line.strip!
          parts = line.split("\t")
          parts[mutation_pos] = translations[parts[mutation_pos]]
          f.puts parts * "\t"
        end
      end
    rescue Exception
      FileUtils.rm filename if File.exists? filename
      raise $!
    ensure
      file.close
    end

    nil
  end


  COSMIC.claim COSMIC.gene_damage_analysis, :proc do

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID",
                    :fields => ["Avg. damage score", "Bg. Avg. damage score", "T-test p-value"],
                    :type => :list, :cast => :to_f, :namespace => COSMIC.organism, :unnamed => true)

    Workflow.require_workflow 'DbNSFP'
    require 'rbbt/util/R'
    database = DbNSFP.database
    database.unnamed = true
    damage_fields = database.fields.select{|f| f =~ /rankscore/ }
    db = COSMIC.knowledge_base.get_database(:gene_principal_isoform_mutations)

    R.eval "a=1" # To start Rserver for all cpus
    RbbtSemaphore.with_semaphore 1 do |sem|
      TSV.traverse db, :cpus => 10, :into => tsv, :bar => "Damage analysis using DbNSFP" do |gene, mis|
        mis = mis.flatten.compact
        next if mis.empty?
        protein = mis.first.partition(":").first
        next unless protein =~ /^ENSP/

        dbNSFP_tsv = database.get_prefix(protein).slice(damage_fields)
        dbNSFP_tsv.unnamed = true

        all_damage_scores = dbNSFP_tsv.collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact
        damage_scores = dbNSFP_tsv.select(:key => mis).collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact

        if damage_scores.length < 3 or all_damage_scores.uniq.length < 3
          damage_score_pvalue = 1
        else
          RbbtSemaphore.synchronize(sem) do
            damage_score_pvalue = R.eval("t.test(#{R.ruby2R(damage_scores)}, #{R.ruby2R(all_damage_scores)}, 'greater')['p.value']").to_f
          end
        end
        [gene, [Misc.mean(damage_scores), Misc.mean(all_damage_scores), damage_score_pvalue]]
      end
    end

    tsv.to_s
  end
end

Misc.add_libdir './lib'
require 'rbbt/sources/COSMIC/indices'
require 'rbbt/entity/COSMIC'
require 'rbbt/knowledge_base/COSMIC'

require 'rbbt'
require 'rbbt/tsv'
require 'rbbt/resource'
require 'rbbt/util/misc/omics'
require 'rbbt/sources/organism'

module COSMIC
  extend Resource
  self.subdir = "share/databases/COSMIC"

  VERSION='v95'
  def self.organism(build)
    Organism.organism_for_build build
  end

  def self.token
    @token ||= Rbbt::Config.get "cosmic_token", "key:cosmic_token", :default => ENV["COSMIC_TOKEN"]
    if @token.nil?
      raise "Please setup you cosmic_token (e.g. in .rbbt/etc/config) following instructions: $(echo 'email@example.com:mycosmicpassword' | base64)"
    else
      @token
    end
  end

  def self.get_real_url(url)
    info = JSON.parse(CMD.cmd("curl -H \"Authorization: Basic #{token}\" #{ url }").read)
    info["url"]
  end

  def self.get_file_url(file, build, version=VERSION)
    url = File.join("https://cancer.sanger.ac.uk/cosmic/file_download/#{build}/cosmic/#{version}/", file)
    get_real_url(url)
  end

  %w(GRCh37 GRCh38).each do |build|
    %w(CosmicMutantExport CosmicStructExport CosmicBreakpointsExport CosmicResistanceMutations CosmicCompleteCNA CosmicCompleteGeneExpression).each do |file|
      COSMIC.claim COSMIC[build][".source"][file + '.gz'], :proc do |filename|
        real_url = COSMIC.get_file_url(file + '.tsv.gz', build)
        Open.mkdir File.dirname(filename)
        CMD.cmd_log("wget -O '#{filename}' '#{real_url}'")
        nil
      end
    end

    #COSMIC.claim COSMIC[build][".source"].cancer_gene_census, :proc do |filename|
    #  url = "sftp-cancer.sanger.ac.uk/cosmic/grch37/cosmic/v82/cancer_gene_census.csv"
    #  raise "Follow #{ url } and place the file uncompressed in #{filename}"
    #end
    #
    #COSMIC.claim COSMIC[build].gene_census, :proc do |filename|
    #  tsv = COSMIC[build][".source"].cancer_gene_census.tsv :key_field => "Entrez GeneId", :fields => [], :type => :single, :sep => ",", :header_hash => ''
    #  res = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => [], :type => :single)
    #  organism = "Hsa/feb2014"
    #  entrez2ensg = Organism.identifiers(organism).tsv :key_field => "Entrez Gene ID", :fields => ["Ensembl Gene ID"], :type => :single
    #  TSV.traverse tsv, :into => res, :bar => true do |entrez|
    #    entrez2ensg[entrez]
    #  end
    #  res.to_s
    #end


    COSMIC.claim COSMIC[build].sample_info, :proc do |file|
      parser = TSV::Parser.new COSMIC[build].mutations, :key_field => "Sample name", :type => :list
      fields = parser.fields
      site_and_histology_fieds = fields.select{|f| f =~ /site|histology/i }
      dumper = TSV::Dumper.new parser.options.merge(:fields => site_and_histology_fieds)
      dumper.init
      pos = site_and_histology_fieds.collect{|f| fields.index f}
      TSV.traverse parser, :bar => true, :into => dumper do |sample, fields|
        [sample, fields.values_at(*pos)]
      end

      TSV.collapse_stream dumper
    end

    COSMIC.claim COSMIC[build].mutations, :proc do |directory|
      Workflow.require_workflow "Sequence"
      url = COSMIC[build][".source"].CosmicMutantExport.produce.find
      job = Sequence.job(:expanded_vcf, "COSMIC", :vcf_file => url)

      stream = Open.open(url)
      parser = TSV::Parser.new stream, :header_hash => "", :key_field => "GENOMIC_MUTATION_ID", :type => :list

      dumper = TSV::Dumper.new parser.options.merge({:fields => parser.fields + ["Genomic Mutation"]})
      dumper.init

      pos_i = parser.identify_field "Mutation genome position"
      cds_i = parser.identify_field "Mutation CDS"
      sample_i = parser.identify_field "ID_sample"
      dumper = TSV.traverse parser, :type => :list, :into => dumper, :bar => url do |mid, values|
        position = values[pos_i]
        cds = values[cds_i]
        sample = values[sample_i]

        mid = [mid, sample] * ":"

        if position.nil? or position.empty?
          Log.debug "Empty genomic position: #{mid}"
          next
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

    COSMIC.claim COSMIC[build].mutation_list, :proc do 
      TSV.traverse TSV.stream_column(COSMIC[build].mutations, "Genomic Mutation"), :type => :array, :into => :stream do |line|
        next if line.split(":").last.include? "?"
        line
      end
    end

    COSMIC.claim COSMIC[build].mi_drug_resistance, :proc do |filename|
      tsv = COSMIC[build][".source"].CosmicResistanceMutations.tsv :key_field => "Transcript", :fields => ["AA Mutation", "Drug Name","Sample ID", "Pubmed Id", "Zygosity"], :type => :double, :merge => true, :header_hash => ''
      res = TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Drug name","Sample","PMID", "Zygosity"], :type => :double)
      organism = "Hsa/feb2014"
      enst2ensp = Organism.transcripts(organism).tsv :key_field => "Ensembl Transcript ID", :fields => ["Ensembl Protein ID"], :type => :single, :merge => true, :persist => true
      TSV.traverse tsv, :into => res, :bar => true do |transcript, values|
        protein = enst2ensp[transcript]
        new = []
        Misc.zip_fields(values).each do |aa_mutation,drug,sample,pmid,zygosity|
          mutation = Misc.translate_prot_mutation_hgvs2rbbt(aa_mutation)
          next if mutation.nil?
          next if protein.nil? or protein.empty?
          mi = [protein,mutation] * ":"
          new << [mi, [drug,sample,pmid,zygosity]]
        end
        new.extend MultipleResult
        new
      end
      res.to_s
    end

    def self.gene2enst(gene, build)
      organism ||= COSMIC.organism(build)
      @@gene2ensg ||= Organism.identifiers(organism).index :target => "Ensembl Gene ID", :persist => true
      @@ensg2enst ||= Organism.transcripts(organism).tsv :key_field => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :type => :flat
      transcripts = case
                  when gene =~ /_ENST/
                    [gene.split("_").last]
                  when gene =~ /ENSG/
                    @@ensg2enst[gene]
                  else
                    ensg = @@gene2ensg[gene] || @@gene2ensg[gene.sub('_HUMAN','')]
                    if ensg.nil?
                      Log.debug "Unknown gene name: #{gene}"
                      nil
                    else
                      @@ensg2enst[ensg]
                    end
                  end
    end

    COSMIC.claim COSMIC[build].completeCNA, :proc do |filename|
      res = TSV::Dumper.new(:key_field => "Sample name", :fields => ["Ensembl Transcript ID","CNA"], :type => :double)
      res.init
      TSV.traverse COSMIC[build][".source"].CosmicCompleteCNA, :key_field => "SAMPLE_NAME", :fields => ["gene_name", "MUT_TYPE"],
        :header_hash => "", :type => :double, :into => res, :bar => true do |sample, values|
        sample = sample.first if Array === sample
        new = []
        Misc.zip_fields(values).each do |gene, cna|
          transcripts = COSMIC.gene2enst(gene)
          next if transcripts.nil?
          cnas = [cna] * transcripts.length
          new << [sample, [transcripts, cnas]]
        end
        new.extend MultipleResult
        new
      end

      TSV.collapse_stream res.stream
    end

    COSMIC.claim COSMIC[build].geneExpression, :proc do |filename|
      res = TSV::Dumper.new(:key_field => "Sample name", :fields => ["Associated Gene Name","Regulation"], :type => :double)
      res.init
      TSV.traverse COSMIC[build][".source"].CosmicCompleteGeneExpression, :key_field => "SAMPLE_NAME", :fields => ["GENE_NAME", "REGULATION"],
        :header_hash => "", :type => :double, :into => res, :bar => true do |sample, values|
        sample = sample.first if Array === sample
        new = []
        Misc.zip_fields(values).each do |gene, expr|
          next if transcripts.nil?
          next if expr == 'normal'
          new << [sample, [gene, expr]]
        end
        new.extend MultipleResult
        new
      end

      TSV.collapse_stream res.stream
    end

    COSMIC.claim COSMIC[build].breakpoints, :proc do |filename|
      res = TSV::Dumper.new(:key_field => "Sample name", :fields => ["Chrom From", "Location From min", "Location From max", "Strand From", "Chrom To", "Location To min", "Location To max", "Strand To"], :type => :double, :namespace => COSMIC.organism(build))
      res.init
      TSV.traverse COSMIC[build][".source"].CosmicBreakpointsExport, :fields => ["Chrom From", "Location From min", "Location From max", "Strand From", "Chrom To", "Location To min", "Location To max", "Strand To"],
        :header_hash => "", :type => :double, :into => res, :bar => true do |sample, values|
        sample = sample.first if Array === sample
        [sample, values]
      end

      TSV.collapse_stream res.stream
    end

    COSMIC.claim COSMIC[build].mutations_hg18, :proc do |filename|
      require 'rbbt/sources/organism'
      organism = COSMIC.organism(build)
      file = COSMIC[build].mutations.open
      begin

        while (line = file.gets) !~ /Genomic Mutation/; end
        fields = line[1..-2].split("\t")
        mutation_pos = fields.index "Genomic Mutation"

        mutations = CMD.cmd("grep -v '^#'|cut -f #{mutation_pos + 1}|sort -u", :in => COSMIC[build].mutations.open).read.split("\n").select{|m| m.include? ":" }

        translations = Misc.process_to_hash(mutations){|mutations| Organism.liftOver(mutations, self.organism(build), "Hsa/may2009")}

        File.open(filename, 'w') do |f|
          f.puts "#: :type=:list#:namespace=#{organism}"
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
        FileUtils.rm filename if File.exist? filename
        raise $!
      ensure
        file.close
      end

      nil
    end

    COSMIC.claim COSMIC[build].gene_damage_analysis, :proc do
      require 'rbbt/knowledge_base/COSMIC'

      tsv = TSV.setup({}, :key_field => "Ensembl Gene ID",
                      :fields => ["Avg. damage score", "Bg. Avg. damage score", "T-test p-value"],
                      :type => :list, :cast => :to_f, :namespace => COSMIC.organism(build), :unnamed => true)

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
end

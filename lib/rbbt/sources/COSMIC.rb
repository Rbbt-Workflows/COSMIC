require 'rbbt'
require 'rbbt/tsv'
require 'rbbt/resource'

module COSMIC
  extend Resource
  self.subdir = "share/databases/COSMIC"

  def self.organism
    Organism.default_code "Hsa"
  end

  COSMIC.claim COSMIC.mutations_register_data, :proc do |filename|
    url = "http://cancer.sanger.ac.uk/files/cosmic/current_release/CosmicMutantExportIncFus.tsv.gz"
    raise "Follow #{ url } and place the file uncompressed in #{filename}"
  end

  COSMIC.claim COSMIC.sample_info, :proc do |file|
    site_and_histology_fieds = TSV.parse_header(COSMIC.mutations).fields.select{|f| f =~ /site|histology/i }
    tsv = COSMIC.mutations.tsv(:key_field => "Sample name", :fields => site_and_histology_fieds, :type => :list)
    tsv.to_s
  end

  COSMIC.claim COSMIC.mutations, :proc do |directory|
    url = COSMIC.mutations_register_data.produce.find
    stream = CMD.cmd('awk \'BEGIN{FS="\t"} { if ($8 != "" && $8 != "Mutation ID") { sub($8, "COSM" $8 ":" $3)}; print}\'', :in => Open.open(url), :pipe => true)

    parser = TSV::Parser.new stream, :header_hash => "", :key_field => "Mutation ID", :type => :list

    dumper = TSV::Dumper.new parser.options.merge({:fields => parser.fields + ["Genomic Mutation"]})
    dumper.init

    pos_i = parser.identify_field "Mutation GRCh37 genome position"
    cds_i = parser.identify_field "Mutation CDS"
    dumper = TSV.traverse parser, :type => :list, :into => dumper, :bar => url do |mid, values|
      position = values[pos_i]
      cds = values[cds_i]

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
          change = case
                   when cds =~ />/
                     cds.split(">").last
                   when cds =~ /del/
                     deletion = cds.split("del").last
                     case
                     when deletion =~ /^\d+$/
                       "-" * deletion.to_i 
                     when deletion =~ /^[ACTG]+$/i
                       "-" * deletion.length
                     else
                       Log.debug "Unknown deletion: #{ deletion }"
                       deletion
                     end
                   when cds =~ /ins/
                     insertion = cds.split("ins").last
                     case
                     when insertion =~ /^\d+$/
                       "+" + "N" * insertion.to_i 
                     when insertion =~ /^[NACTG]+$/i
                       "+" + insertion
                     else
                       Log.debug "Unknown insertion: #{insertion }"
                       insertion
                     end
                   else
                     Log.debug "Unknown change: #{cds}"
                     "?(" << cds << ")"
                   end
          [mid, values + [position + ":" + change]]
        end
      end
    end
    dumper.stream
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
        [gene, [Misc.mean(all_damage_scores), Misc.mean(damage_scores), damage_score_pvalue]]
      end
    end

    tsv.to_s
  end
end

Misc.add_libdir './lib'
require 'rbbt/sources/COSMIC/indices'
require 'rbbt/sources/COSMIC/entity'
require 'rbbt/sources/COSMIC/knowledge_base'

if __FILE__ == $0
  require 'rbbt/workflow'
  ppp COSMIC.sample_info.produce(true).tsv.to_s
end

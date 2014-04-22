require 'rbbt'
require 'rbbt/tsv'
require 'rbbt/resource'

module COSMIC
  extend Resource
  self.subdir = "share/databases/COSMIC"

  COSMIC.claim COSMIC.mutations, :proc do 
    url = "ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/CosmicCompleteExport_v68.tsv.gz"

    stream = CMD.cmd('awk \'BEGIN{FS="\t"} { if ($12 != "" && $12 != "Mutation ID") { sub($12, "COSM" $12 ":" $4)}; print}\'', :in => Open.open(url), :pipe => true)
    tsv = TSV.open(stream, :type => :list, :header_hash => "", :key_field => "Mutation ID", :namespace => "Hsa/jun2011")
    tsv.fields = tsv.fields.collect{|f| f == "Gene name" ? "Associated Gene Name" : f}
    tsv.add_field "Genomic Mutation" do |mid, values|
      position = values["Mutation GRCh37 genome position"]
      cds = values["Mutation CDS"]

      if position.nil? or position.empty?
        nil
      else
        position = position.split("-").first

        chr, pos = position.split(":")
        chr = "X" if chr == "23"
        chr = "Y" if chr == "24"
        chr = "M" if chr == "25"
        position = [chr, pos ] * ":"

        if cds.nil?
          position
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
          position + ":" + change
        end
      end
    end

    tsv.to_s.gsub(/(\d)-(\d)/,'\1:\2')
  end

  COSMIC.claim COSMIC.mutations_hg18, :proc do |filename|
    require 'rbbt/sources/organism'
    file = COSMIC.mutations.open
    begin

      while (line = file.gets) !~ /Genomic Mutation/; end
      fields = line[1..-2].split("\t")
      mutation_pos = fields.index "Genomic Mutation"

      mutations = CMD.cmd("grep -v '^#'|cut -f #{mutation_pos + 1}|sort -u", :in => COSMIC.mutations.open).read.split("\n").select{|m| m.include? ":" }

      translations = Misc.process_to_hash(mutations){|mutations| Organism.liftOver(mutations, "Hsa/jun2011", "Hsa/may2009")}

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
    Workflow.require_workflow "MutEval"
    require 'db_nsfp'
    require 'rsruby'

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", 
                    :fields => ["Avg. damage score", "Bg. Avg. damage score", "T-test p-value"], 
                    :type => :list, :cast => :to_f, :namespace => COSMIC.organism, :unnamed => true)

    database = DbNSFP.database
    database.unnamed = true
    damage_fields = database.fields.select{|f| f =~ /converted/ }
    db = COSMIC.knowledge_base.get_database(:gene_principal_isoform_mutations)
    db.unnamed = true
    db.with_monitor :desc => "Damage analysis using DbNSFP", :step => 10000 do
      db.through do |gene,mis|
        next if mis.empty?
        protein = mis.first.partition(":").first
        next unless protein.index "ENSP"

        dbNSFP_tsv = database.get_prefix(protein).slice(damage_fields)
        dbNSFP_tsv.unnamed = true

        all_damage_scores = dbNSFP_tsv.collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact
        damage_scores = dbNSFP_tsv.select(:key => mis).collect{|k,values| good = values.reject{|v| v == -999}; good.any? ? Misc.mean(good) : nil}.compact

        if damage_scores.length < 3 or all_damage_scores.uniq.length < 3
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

require 'rbbt/sources/COSMIC/indices'
require 'rbbt/sources/COSMIC/entity'
require 'rbbt/sources/COSMIC/knowledge_base'


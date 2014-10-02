Log.deprecated "do not require #{__FILE__}, require 'rbbt/entity/COSMIC' instead. In: #{caller.select{|l| l =~ /rbbt|workflow/}.first}"

require 'rbbt/entity/COSMIC'

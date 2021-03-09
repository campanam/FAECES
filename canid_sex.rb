#!/usr/bin/env ruby
require 'bigdecimal'
require 'ostruct'
require 'optparse'
#----------------------------------------------------------------------------
# Calculate binomial probability
def binom(n, k, p = 0.50)
	if n - k > k # Maximize efficient division
		nCk = factorial(n, n-k)
		nCk /= factorial(k)
	else
		nCk = factorial(n, k) # Calculate factorial(n)/factorial(k)
		nCk /= factorial(n - k)
	end
	prob = BigDecimal(p.to_s) ** BigDecimal(k.to_s) * BigDecimal((1.0 - p).to_s) ** BigDecimal((n-k).to_s) * BigDecimal(nCk.to_s)
	return prob
end
#----------------------------------------------------------------------------
# Calculate factorials
def factorial(n, k = 1)
	total = 1
	while n > k # Permits efficient division of factorials
		total *= n
		n -= 1
	end
	return total
end
#----------------------------------------------------------------------------
def call_sex(falsepos, falseneg, n, y)
	if n == 0 or n < $options.minalleles # Insurance for invalid values of minalleles like -1, 0, etc.
		sex = "?"
	elsif y == 0 && falseneg > $options.alpha
		sex = "?"
	elsif y == 0 && falseneg <= $options.alpha
		sex = "F"
	elsif y > 0 && falsepos > $options.alpha
		sex = "M"
	elsif y > 0 && falsepos <= $options.alpha
		sex = "F"
	end
	return sex
end
#----------------------------------------------------------------------------
def build_species_hash
	$species_hash = {} # Hash assigning individuals to species. Not defined if file not provided
	File.open($options.species) do |f3|
		while line = f3.gets
			sam_arr = line.strip.split
			$species_hash[sam_arr[0]] = $species_hash[sam_arr[1]]
		end
	end
end
#----------------------------------------------------------------------------
def build_Y_hash
	# Hashes of alleles in form in [X allele, Y allele]
	$kf_alleles = { '19789242' => ['C','T'], '19789276' => ['T','G'], '19789374' => ['C','G'], '19789494' => ['C','G'] } # Y alleles only in kit fox
	$sex_alleles = { '19789218' => ['G','A'], '19789227' => ['A','G'], '19789233' => ['G','A'],
				'19789237' => ['G','T'], '19789239' => ['C','T'], '19789251' => ['C','G'],
				'19789275' => ['T','C'], '19789293' => ['G','A'], '19789314' => ['G','A'],
				'19789317' => ['A','G'], '19789321' => ['G','A'], '19789356' => ['C','T'],
				'19789401' => ['C','T'], '19789437' => ['C','T'], '19789440' => ['C','T'],
				'19789458' => ['A','G'], '19789503' => ['G','A'], '19789515' => ['G','A'],
				'19789551' => ['G','A'], '19789563' => ['A','G'], '19789569' => ['T','C'] }
end
#----------------------------------------------------------------------------
def calculate_Y_frequency
	y_sum = 0.0 # Sum of allele frequencies to calculate mean/Sum of Y alleles
	denom = 0.0 # Number of known males with estimated values/Number of total alleles
	ymin = 1.0 # Minimum Y frequency
	stdev = [] # Array of vals to calculate stdev for z-score
	stdevxyz = [] # Array of vals to calculate stdev for z-score of Y-X ratio
	File.open($options.males) do |f2|
		while line = f2.gets
			male = line.strip
			unless $sample_sex[male].nil?
				if $sample_sex[male].sum >= $options.minalleles  && $sample_sex[male][1] > 0 # Filter out low quality samples/false negatives
					if $options.mean
						denom += 1.0
						yval = $sample_sex[male][1].to_f/$sample_sex[male].sum.to_f
						y_sum += yval
						unless $options.zscore.nil?
							stdev.push(yval)
						end
						stdevxyz.push($sample_sex[male][1].to_f/$sample_sex[male][0].to_f) if $options.xyz != ""
					elsif $options.min
						yval = $sample_sex[male][1].to_f/$sample_sex[male].sum.to_f
						ymin = yval if yval < ymin
					else # Instead of using mean sample frequency, calculate overall rate based on total allele counts
						y_sum += $sample_sex[male][1].to_f
						denom += $sample_sex[male].sum.to_f
					end
				end
			end
		end
	end
	if $options.min 
		$options.yfreq = ymin
		$stderr.puts "Empirical minimum Y-frequency: " + ymin.to_s
	else
		$options.yfreq = y_sum/denom
	end
	$stderr.puts "Empirical mean Y-frequency: " + (y_sum/denom).to_s
	xyratio = (y_sum/(denom - y_sum))
	$stderr.puts "Empirical mean Y-X ratio: " + xyratio.to_s
	unless $options.zscore.nil?
		total = 0.0
		for sd in stdev
			total += (sd - $options.yfreq) ** 2
		end
		sd = Math.sqrt(total / (stdev.size - 1))
		$options.yfreq += $options.zscore * sd
		$stderr.puts "Empirical Y-frequency StDev: " + sd.to_s
	end
	if $options.xyz != ""
		total = 0.0
		for sd in stdevxyz
			total += (sd - xyratio) ** 2
		end
		sd = Math.sqrt(total / (stdevxyz.size - 1))
		xyratio += $options.xyz * sd
		$stderr.puts "Empirical Y-X ratio StDev: " + sd.to_s
		$options.xyratio = xyratio
		$stderr.puts "Empirical Y-X ratio cut-off: " + $options.xyratio.to_s
	else
		$stderr.puts "Empirical Y-frequency cut-off: " + $options.yfreq.to_s
	end
end
#----------------------------------------------------------------------------
def process_vcf
	start = false
	File.open($options.infile) do |f1|
		while line = f1.gets
			if line[0..5] == "#CHROM"
				start = true
				@samples = line.strip.split("\t")[9..-1]
				$sample_sex = {} # Hash of X chromosome allele, Y alleles counts
				for sample in @samples
					$sample_sex[sample] = [0,0]
				end
			elsif start
				line_arr = line.strip.split("\t")
				if line_arr[0] == "CM000039" or line_arr[0] == "Coy_ZF_X" or line_arr[0] == "Coy_ZF_Y" or line_arr[0] == "KF_ZF_X" or line_arr[0] == "KF_ZF_Y"
					sexsite = line_arr[1]
					sexsite = (sexsite.to_i + 19789197).to_s unless line_arr[0] == "CM000039"
					if $sex_alleles.keys.include?(sexsite)
						alleles = [line_arr[3], line_arr[4].split(",")].flatten
						# Get AD site
						format = line_arr[8].split(":")
						ad_index = format.index("AD")
						x_index = alleles.index($sex_alleles[sexsite][0])
						y_index = alleles.index($sex_alleles[sexsite][1])
						# If reference X or Y alleles are not identified in dataset, will return nil.
						for sample in @samples
							sample_alleles = line_arr[@samples.index(sample) + 9].split(":")[ad_index].split(",")
							$sample_sex[sample][0] += sample_alleles[x_index].to_i unless x_index.nil?
							$sample_sex[sample][1] += sample_alleles[y_index].to_i unless y_index.nil?
						end
					elsif $kf_alleles.keys.include?(sexsite) and !$species_hash.nil?
						alleles = [line_arr[3], line_arr[4].split(",")].flatten
						# Get AD site
						format = line_arr[8].split(":")
						ad_index = format.index("AD")
						x_index = alleles.index($kf_alleles[sexsite][0])
						y_index = alleles.index($kf_alleles[sexsite][1])
						# If reference X or Y alleles are not identified in dataset, will return nil.
						for sample in @samples
							if $sample_hash[sample] == "KF" # Some ugly hard-coding for the ages.
								sample_alleles = line_arr[@samples.index(sample) + 9].split(":")[ad_index].split(",")
								$sample_sex[sample][0] += sample_alleles[x_index].to_i unless x_index.nil?
								$sample_sex[sample][1] += sample_alleles[y_index].to_i unless y_index.nil?
							end
						end
					end
				end
			end
		end
	end
	calculate_Y_frequency if $options.males != ""
	if $options.xyratio != ""
		puts "Sample\tX-alleles\tY-alleles\tY-XRatio\tCalledSex"
	else
		puts "Sample\tX-alleles\tY-alleles\tFalsePosY\tFalseNegY\tCalledSex"
	end
	for sample in @samples
		n = $sample_sex[sample].sum
		if n > 0
			if $sample_sex[sample][1] == 0
				falsepos = "NA"
				if $sample_sex[sample][1] == 0
					falseneg = binom(n, 0, $options.yfreq)
				end
				falseneg = falseneg.round(3)
			else
				falseneg = "NA"
				# We are testing for Y drop-in rather than X drop-out
				falsepos = 0.0
				for k in 1 .. $sample_sex[sample][1]
					falsepos += binom(n, k, $options.yfreq)
				end
				falsepos = falsepos.round(3)
			end
		else
			falsepos = falseneg = "NA"
		end
		if $options.xyratio != ""
			xyratio = $sample_sex[sample][1].to_f/$sample_sex[sample][0].to_f
			if n < $options.minalleles
				sex = "?"
			elsif xyratio < $options.xyratio
				sex = "F"
			else
				sex = "M"
			end
			puts sample + "\t" + $sample_sex[sample].join("\t") + "\t" + xyratio.to_s + "\t" + sex
		else
			sex = call_sex(falsepos, falseneg,n, $sample_sex[sample][1])
			puts sample + "\t" + $sample_sex[sample].join("\t") + "\t" + falsepos.to_s + "\t" + falseneg.to_s + "\t" + sex
		end
	end
end
#----------------------------------------------------------------------------
class Parser
	def self.parse(options)
		args = OpenStruct.new
		args.infile = "" #Input VCF
		args.species = "" # File of species assignments
		args.males = "" # List of known males
		args.mean = false # Infer male Y frequency using mean value across samples
		args.zscore = nil # Infer male Y frequency using z-scores
		args.min = false # Infer minimum Y frequency using minimum value across samples
		args.minalleles = 10 # Minimum number of allele reads to call sex
		args.yfreq = 0.5 # Y allele frequency
		args.alpha = 0.05 # Alpha value
		args.xyratio = "" # Use raw X-Y allele ratio to estimate sex
		args.xyz = "" # Use z-score cut-off to determine X-Y allele ratio
		opt_parser = OptionParser.new do |opts|
			opts.banner = "Command-line usage: ruby canid_sex.rb [options]"
			opts.separator ""
			opts.on("-i","--input [FILE]", String, "Input VCF file") do |vcf|
				args.infile = vcf
			end
			opts.on("-s","--species [FILE]", String, "Species assignment file") do |vcf|
				args.species = vcf if vcf != nil
			end
			opts.on("-k", "--males [FILE]", String, "List of known males to empirically calculate Y-allele frequency") do |males|
				args.males = males if males != nil
			end
			opts.on("--mean", "Infer Y frequency using mean frequency across samples") do
				args.mean = true
				args.min = false
			end
			opts.on("-z", "--zscore [VALUE]", Float, "Infer mean Y frequency and set z-score bound") do |zsc|
				if zsc != nil
					args.mean = true
					args.min = false
					args.zscore = zsc
				end
			end
			opts.on("--min", "Infer Y frequency using minimum frequency across samples") do
				args.min = true
				args.mean = false
			end
			opts.on("-m", "--minalleles [VALUE]", Integer, "Minimum number of allele reads to call sex (Default = 10)") do |sex|
				args.minalleles = sex if sex != nil
			end
			opts.on("-y", "--yfreq [VALUE]", Float, "Expected Y-allele frequency (Default = 0.5)") do |alph|
				args.yfreq = alph if alph != nil
			end
			opts.on("-a", "--alpha [VALUE]", Float, "Alpha value (Default = 0.05)") do |alph|
				args.alpha = alph if alph != nil
			end
			opts.on("--yx [VALUE]", Float, "Use specified Y-X allele ratio to determine sex rather than statistical test.") do |xyratio|
				args.xyratio = xyratio if xyratio != nil
			end
			opts.on("-Z [VALUE]", Float, "Use z-score cut-off to determine Y-X ratio.") do |xyz|
				args.xyz = xyz if xyz != nil
				args.zscore = xyz
				args.mean = true
				args.min = false # Set this all up, args.xyratio turned on later once calculated in calculate_Y_frequency
			end
			opts.on_tail("-h","--help", "Show help") do
				puts opts
				exit
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#----------------------------------------------------------------------------
BigDecimal.limit(16)
ARGV[0] ||= "-h"
$options = Parser.parse(ARGV)
build_Y_hash
if !FileTest.exist?($options.infile)
	puts "Input file not found. Exiting."
else
	if !FileTest.exist?($options.species)
		$stderr.puts "Species assignment file not found. Ignoring species assignments."
	else
		build_species_hash
	end
	process_vcf
end
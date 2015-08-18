#!/usr/bin/env ruby

def parse_fname(fname)
  { dir: File.dirname(fname),
    base: File.basename(fname, File.extname(fname)),
    ext: File.extname(fname) }
end

def get_end_posn(read_len, start_posn)
  start_posn + read_len -1
end

require 'parse_fasta'
require 'trollop'

Signal.trap("PIPE", "EXIT")

opts = Trollop::options do
  banner <<-EOS

  Lalalala!

  Make perfectly perfect read pairs from multiple genomes.

  --coverage-file
    name\tcov

  Lalalala!

  Options:
  EOS
  opt :fasta, 'Fasta file name', type: :string
  opt(:coverage_file, 'The file with the coverages', type: :string)
  opt :read_len, 'The length of reads.', type: :int, default: 150
  opt :insert_len, 'Insert length', type: :int, default: 100
  opt :qual_score, 'The fake qual score', type: :string, default: "I"
  opt(:outdir, 'Output directory', type: :string,
      default: '/home/cshatley/projects/viral_abundance/reads')
end

if opts[:fasta].nil?
  Trollop.die :fasta, "You must enter a file name"
elsif !File.exist? opts[:fasta]
  Trollop.die :fasta, "The file must exist"
end

if opts[:coverage_file].nil?
  Trollop.die :coverage_file, "You must enter a file name"
elsif !File.exist? opts[:coverage_file]
  Trollop.die :coverage_file, "The file must exist"
end

# Given a pre_read and an insert_len, split the pre-read into 2 actual
# reads, returns
def split pre_read, read_len
  [pre_read[0, read_len],
   pre_read.reverse[0, read_len].tr('ACTG', 'TGAC')]
end

fname = parse_fname opts[:fasta]
f1 = File.join(opts[:outdir],
               "#{fname[:base]}." +
               "len_#{opts[:read_len]}.ins_#{opts[:insert_len]}.1.fq")
f2 = File.join(opts[:outdir],
               "#{fname[:base]}." +
               "len_#{opts[:read_len]}.ins_#{opts[:insert_len]}.2.fq")
out1 = File.open f1, "w"
out2 = File.open f2, "w"

coverages = {}
File.open(opts[:coverage_file]).each_line do |line|
  header, cov = line.chomp.split "\t"

  if coverages.has_key? header
    abort "#{header} is repeated in #{opts[:coverage_file]}"
  else
    coverages[header] = cov.to_i
  end
end

read_len = opts[:read_len] * 2 + opts[:insert_len]
orig_len = opts[:read_len] # for use in splitting an coverage

frag_num = 1
FastaFile.open(opts[:fasta]).each_record do |header, sequence|
  unless coverages.has_key? header
    abort("#{header} is present in #{opts[:fasta]} but not in " +
          "#{opts[:coverages_file]}")
  end

  the_cov = coverages[header]
  cov = the_cov * (1 / (orig_len*2 / read_len.to_f))

  the_step = (read_len / cov.to_f).round
  start_posns = (0..read_len).step(the_step)
  start_posns =
    start_posns.take(cov)

  genome_len = sequence.length

  start_posns.each do |posn|
    number_to_take =
      ((genome_len - posn + 1) / read_len.to_f).floor

    these_start_posns =
      (posn..genome_len).step(read_len).take(number_to_take).to_a

    these_start_posns.each do |start_posn|
      end_posn = get_end_posn(read_len, start_posn)
      this_fragment = sequence[start_posn..end_posn]

      read1, read2 = split this_fragment, orig_len

      out1.printf "@#{header}_read_%s/1\n", frag_num
      out1.puts read1
      out1.printf("+%s frag_start=%d frag_end=%d\n",
             header, start_posn+1, end_posn+1)
      out1.puts opts[:qual_score] * read1.length

      out2.printf "@#{header}_read_%s/2\n", frag_num
      out2.puts read2
      out2.printf("+%s frag_start=%d frag_end=%d\n",
             header, start_posn+1, end_posn+1)
      out2.puts opts[:qual_score] * read2.length

      frag_num += 1
    end
  end
end

out1.close
out2.close

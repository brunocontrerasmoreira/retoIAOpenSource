#!/usr/bin/env perl
use strict;
use warnings;

# Bruno Contreras Nov2025

if(!defined($ARGV[1])) {
  die "# usage: $0 seqs.fa classname\n";
}  
 
my ($fastafile,$class) = @ARGV;

my ($id,%seq);
open(FA,"<",$fastafile) || 
  die "# ERROR: cannot read $fastafile\n";
while(<FA>) {
  if(/^>(\S+)/) {
    $id = $1;

  } else {
    chomp;
    $seq{ $id } .= uc($_);
  }
}
close(FA);

print "id;sequence;ChIPseq;AA;AC;AG;AT;CA;CC;CG;GA;GC;TA\n";
foreach $id (sort keys(%seq)) {

  print "$id;$seq{$id};$class";
  print dinucleotide_frequencies($seq{$id}) . "\n";
}


# modified from Claude code in response to queries:
# i) can you produce a Perl function that takes a string of nucleotides and return the frequencies of dinucleotides?
# ii) can you modify it so that reverse dinucleotides are merged with forward ones?

sub reverse_complement_dinuc {
    my ($dinuc) = @_;
    my %complement = ('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G');
    my $rev = reverse($dinuc);
    my $rc = join('', map { $complement{$_} } split('', $rev));
    return $rc;
}

sub canonical_dinuc {
    my ($dinuc) = @_;
    my $rc = reverse_complement_dinuc($dinuc);
    # Return the lexicographically smaller one as canonical
    return ($dinuc lt $rc) ? $dinuc : $rc;
}

# returns CSV frequency string
sub dinucleotide_frequencies {
    my ($sequence) = @_;

    my ($nt1,$nt2,$dinuc,$canonical,%dinuc_count);
    my $freqs = '';
    my $len = length($sequence);

    if($len < 2) {
        return "0;0;0;0;0;0;0;0;0;0"
    }    

    # Count dinucleotides with canonical form
    for (my $i = 0; $i < $len - 1; $i++) {
        $dinuc = substr($sequence, $i, 2);
        $canonical = canonical_dinuc($dinuc);
        $dinuc_count{$canonical}++;
    }

    # Calculate frequencies (normalize by total dinucleotides)
    my $total = $len - 1;
    my $total_freq = 0;
    my %dinuc_freq;

    foreach $nt1 ('A','C','G','T') {
        foreach $nt2 ('A','C','G','T') {
            $dinuc = $nt1.$nt2;
            $canonical = canonical_dinuc($dinuc);

            next if(defined($dinuc_freq{$canonical}));
            #print ";$canonical";

            if(defined($dinuc_count{$canonical})) {	    
                $dinuc_freq{$canonical} = sprintf("%1.3f",$dinuc_count{$canonical} / $total || 0);
	    } else {
                $dinuc_freq{$canonical} = '0';
	    }

	    $freqs .= ";$dinuc_freq{$canonical}";
	    $total_freq += $dinuc_freq{$canonical};  
        }
    } #print "$total_freq\n";

    return $freqs;
}

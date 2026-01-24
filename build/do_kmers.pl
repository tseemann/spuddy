#!/usr/bin/env perl
use 5.32.0;
use strict;
use Path::Tiny;
use List::Util qw(shuffle);
use Data::Dumper;

srand(42); # same 'shuffle' results each time

my $ZSTD = 'zstd -T\$CPUS';
my $PIGZ = 'pigz -p \$CPUS';
my $GZIP = "libdeflate-gzip";
my $ZCAT = "$GZIP -d -c -f";
my $CACHE = "vmtouch -t";

my $GTDB = "/home/shared/db/gtdb/r226/gtdb_genomes_r226";
@ARGV or die "$0 KMER STEP MIN LIM";
my $KMER = shift(@ARGV) // 7;
my $STEP = shift(@ARGV) // 1;
my $MIN = shift(@ARGV) // 10;
my $LIM = shift(@ARGV) // 1000;

say STDERR "Generating: K=$KMER S=$STEP MIN=$MIN LIM=$LIM";
my $name = "kmers.K$KMER.S$STEP.M$MIN.L$LIM";
my $suffix = "K$KMER.S$STEP.gz";
say STDERR "Prefix: $name";
say STDERR "Suffix: $suffix (for re-usable kmers)";

#     0	Representative genome
#     1	GTDB species
#     2	GTDB taxonomy
#     8	No. clustered genomes
#     9	Clustered genomes (CSV)

my %tax_of;
my $taxid = 0;
my @groups;
my @taxa;
my $DIV = 'Bacteria';
foreach (path('sp_clusters.tsv')->lines) {
  chomp;
  my @col = split m/\t/;  
  next unless $col[8] >= $MIN;
  next unless $col[2] =~ m/d__$DIV/;
  my($genus,$species) = split ' ', $col[1];
  $genus =~ s/s__//;
  my $gs = "$DIV/$genus/$species";
  if ($genus =~ m/[0-9]/) {
    say STDERR "Skipping: $gs";
    next;
  }
  path($gs)->mkpath;
  push @groups, "$gs\n";
  ++$taxid;
  $tax_of{$taxid} = $gs;
  push @taxa,"$taxid\t$gs\n";
  say STDERR "$taxid $col[8] $gs";
}
my $G = "$name.groups";
path($G)->spew(\@groups);
path("$name.taxa")->spew(\@taxa);

# main shell script
open my $OUT, '>', "$name.sh";
select $OUT;

#say "set -uex -o pipefail"; ## won't push through
say "#!/usr/bin/env bash";
say "set -eu -o pipefail";
say "CSVTK_THREADS=1";
say "SEQKIT_THREADS=1";
say "LC_ALL=C"; 
say 'CPUS=$(( $(nproc) / 2 ))'; # hope we have > 1

my $LOGFILE = "$name.sh.log";
say "rm -f $LOGFILE";
# The + means append to the log file
my $PARA = "parallel --bar --joblog +$LOGFILE";
my $PARAJ = "$PARA -j \$CPUS";

banner("Making folders");
say "$PARAJ -a $G mkdir -p";
my @G = sort(path($G)->lines({chomp=>1}));

my @faa;
my @gcmd;
#my @kfile;
for my $g (@G) {
  say STDERR "Prepping: $g";
  my @kmerf;
  for my $a (shuffle  path("$GTDB/$g/proteomes.txt")->lines({chomp=>1}) ) {
    push @faa, "$g/$a";
    push @kmerf, "$a.$suffix";
    last if @kmerf >= $LIM; 
  }

  push @gcmd, 
    species_pepmer( $g, $name, scalar(@kmerf) );

  # keep track of subset of files we kmered
  path("$g/$name.kmers.fofn")->spew(
    map { "$_\n" } @kmerf);
}

my $faa_fofn = path("$name.faa.fofn");
$faa_fofn->spew(map { "$_\n" } @faa);

banner("Generating kmers for every isolate");
say "$PARAJ -a $faa_fofn ".dq(clean("
    [[ -f '{}.$suffix' ]] ||
      seqkit sliding -W $KMER -s $STEP $GTDB/{}
    | tr '[:lower:]' '[:upper:]' 
    | grep -v -E '^>|[XUOW*]'
    | tsvtk uniq -H -o {}.$suffix
  "));

# count uniq kmers per isolate
my $gcmd = path("$name.gcmd");
$gcmd->spew(map { "$_\n" } @gcmd);
banner("Finding unique kmers per isolate");
say "$PARAJ -a $gcmd";

# compile all the uniqmers
banner("Identifying global singleton kmers (slow)");
#say "$PARA -j 1 -a $G $CACHE {}/$name";
my $NG = scalar(@G)+1;
#  [[ -r '$name.kmers' ]] ||
say clean("
  xargs -a $G -I % echo %/$name 
  | tr '\\n' '\\0'
  | LC_ALL=C sort --merge --files0-from=-
                  --batch-size=$NG
                  --parallel=\$CPUS
  | LC_ALL=C uniq -u
  > $name.kmers
");

# For each species, find their species-specific ones
say "$CACHE $name.kmers";
banner("Identifying kmers unique to each speices");
#    [[ -r '{}/$name' ]] ||
say "$PARAJ -a $G ".dq(clean("
    LC_ALL=C sort --merge $name.kmers {}/$name
    | LC_ALL=C uniq -d 
    > {}/$name.uniq
    "));

# combine final uniqmers with our porportion data
banner("Add scores to species-specific kmers");
say "$PARAJ -a $G ".dq(clean("
       tsvtk join -H --left-join
             {}/$name.uniq {}/$name.freq
    |  tsvtk mutate2 -H -e \\'{}\\'
             -o {}/$name.anno
    "));

# collate stats of uniqmers for each species
banner("Collect some species stats");
say "$PARAJ -a $G 'wc -l {}/$name.uniq' | sort -bnr > $name.counts";

# collate all the species kmers/socre
banner("Make final databases");
# the -k is important here to ensure it remians sorted?
say clean("
  xargs -a $G -I % echo %/$name.anno
  | tr '\\n' '\\0'
  | LC_ALL=C sort --merge --files0-from=-
                  --batch-size=$NG
                  -k 1.1,1.$KMER
  | tsvtk join -H -f '4;2' - $name.taxa
  | tsvtk cut -H -f1,2,5
  | tsvtk replace -H -f 2 -p '(\\.\\d*[1-9])0+\$|\\.0+\$' -r '\$1'
  > $name.final
");

# generate cut down versions
banner("Shrink the databases");
say "$CACHE $name.final";
say "LINES=\$(wc -l < $name.final)"; # should do in sort command with tee
for my $PC (qw"01 50 75 90 95 99") {
  my $stem = "$name.db.${PC}pc";
  say clean("
    (   cat $name.final 
      | pv -l -s \$LINES -pec -N $stem
      | awk '\$2 >= 0.${PC}0'
      > $stem 
      ) &
  ");
}
say("wait");
sat("rm -fv $name.db.*");
say("zstd --keep $name.db.*"

banner("Cleanup");
say "rm -fv $name.{faa.fofn,gcmd}";

# show the results
banner("Results");
say "wc -l $name.db*pc";
say "\\ls -1shd $name.db*";

# finish up
banner("Pipeline completed!");
say "echo  less -S $name.counts";

system "chmod +x $name.sh";
say STDERR "Run this: ./$name.sh";
exit(0);

#...........................................
sub species_pepmer {
  my($dir,$name,$num) = @_;
  my $outname = "$dir/$name";
#    [[ -f '$outname' ]] ||
  return clean("
    xargs -a $dir/$name.kmers.fofn 
          -I % $ZCAT $dir/%
    | tsvtk freq -H --sort-by-key
    | tsvtk mutate2 -H --at 2
            -w 6 -e '\$2 / $num'
    | tee $dir/$name.freq
    | cut -f 1 
    > $outname
  ");
}
#...........................................
sub dq { 
  $_ = shift;
  s/"/\\"/g;
  return qq/"$_"/;
}
#...........................................
sub clean {
  return join ' ', 
  split ' ', join ' ', @_;
}
#...........................................
sub banner {
  my($text) = @_;
  my $width = $ENV{'COLUMNS'} || 40; 
  my $line = '='x($width-2);
  say "echo '$line'";
  say "echo '$text'";
  #say "echo '$line'";
}
#...........................................

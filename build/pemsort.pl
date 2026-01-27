#!/usr/bin/env perl
use 5.32.0;
use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use File::Spec;
#use File::Which;
use File::Temp qw(tempdir);

my $EXE = basename($0);
my $VER = '0.2.0';

my %opt = (
  'j'=>0, 'k'=>3, 'o'=>'', 'f'=>'', 'S'=>'',
  'T'=>$ENV{'TMPDIR'} // '',
);

sub msg { say STDERR "@_" unless $opt{'q'} };
sub err { msg("ERROR:", @_); exit(1); }

getopts('j:k:f:T:o:vhqsS:', \%opt);

$opt{v} && do { msg("$EXE $VER"); exit(0); };
$opt{h} and show_help(0);
is_int($opt{k}, '-k INT', 2, 1024);
-f $opt{o} and err("Output -o $opt{o} already exists");

my $cpus = cpus();
is_int($opt{j}, '-j CPUS', 0, $cpus);
$opt{j} ||= $cpus; # 0 = all CPUs
my $ram = ram_gb();

my @file = @ARGV;
push @file, read_lines($opt{f}) if $opt{f};
@file >= 2 or err("Need at least 2 input files to sort!");
@file = map { File::Spec->rel2abs($_) } @file;

msg("This is $EXE $VER");
msg("System has $cpus CPUs and $ram GB RAM");
msg("Have",0+@file,"files to sort.");    

my $TD = tempdir(CLEANUP=>1);
msg("Working in: $TD");
my $id=0;

open my $MF, '>', "$TD/Makefile";
select $MF;

say ".PHONY: all";
say "RM := rm -f";
say "SORT := sort $opt{S} --merge -T $TD --batch-size=$opt{k}";

if ($opt{'s'}) {
  # FIXME
  # inputs are not already sorted
  # so we need to sort them first
  msg("-s | pre-sorng not supported yet");
}

if ($opt{'z'}) {
  msg("-z | ompression not supported yet");
}

my %is_temp;

while (@file > 1) {
  #say "# FILES: @file";
  my @in = splice @file, 0, $opt{k};
  $id++;
  my $out = "$TD/pemsort.$id";
  $is_temp{$out}++;
  my @rm = grep { $is_temp{$_} } @in;
  say "$out : @in";
  say "\t\$(SORT) \$^ > \$@";
  say "\t\$(RM) @rm" if @rm;
  push @file, $out;
}
#say "# FILES: @file";
@file==1 or err("Expected only 1 file at end!");
say "all : @file";
if ($opt{o}) {
  say "\tmv -f \$< $opt{o}";
}
else {
  say "\tcat \$< >&3";
  say "\t\$(RM) \$<";
}
close $MF;

my $mopt = $opt{q} ? '--silent' : '';
system('bash', '-c', "make $mopt -j $opt{j} -C $TD all 3>&1 1>&2")==0
  or err("Error sorting files: $!");

msg("$TD should be cleaned up now.");
msg("Sorting compelte.");
exit(0);
#.............................................
sub ram_gb {
  my($line) = grep { m/Mem/ } qx"free -wg";
  $line =~ m/(\d+)/ or err("Can't determine free RAM");
  return $1;
}
#.............................................
sub cpus {
  my($num) = qx"nproc";
  $num or err("Can't determine number of CPUs");
  chomp $num;
  return $num // 1;
}
#.............................................
sub read_lines {
  my($filename) = @_;
  open my $FILE, '<', $filename
    or err("Can't open '$filename'");
  my @lines = <$FILE>;
  chomp @lines;
  return @lines;  
}
#.............................................
sub is_int {
  my($x, $name, $min, $max) = @_;
  !defined($x) 
  || $x !~ m/^\d+$/
  || $x < $min
  || $x > $max
  and err("$name must be ineger betweeb $min and $max");
  return $x;
}
#.............................................
sub show_help {
  my($exitcode) = @_;
  print <<"END_HELP";
SYNOPSIS
  $EXE $VER - parallel external merge sort
USAGE
  $EXE [opts] [-f FOFN | FILES ....]
OPTIONS
  -h       Show this help 
  -v       Show version and exit
  -q       Quiet mode
  -j CPUS  Theeads to use
  -k INT   Files to merge per CPU [$opt{k}]
  -f FOFN  FIle of filenames to sort
  -T DIR   Fast temporary directory [$opt{T}]
  -n       Dry-run - just print steps
  -z       Use LZ4 for intermeidate files
  -s       Sort input files first
  -S STR   Options to pass to Unix 'sort'
END  
END_HELP
  exit($exitcode);
}

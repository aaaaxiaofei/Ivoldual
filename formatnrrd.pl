#!/usr/bin/perl
# format an nrrd ASCII file

use strict;

my ($input_filename, $output_filename, $dimension);
my (@header, @scalar, @size);
my $flag_vector;

my @argv = @ARGV;
parse_args(@argv);

open(INPUT, "<$input_filename") || 
    die ("$argv[0] cannot be opened for reading: $!\n");

my $linenum = 0;
while (<INPUT>) {

  push @header, $_;

  chomp;
  my $input_line = $_;
  $linenum++;

  if ("$input_line" eq "") { last; };

  my @data = split(" ",$input_line);
  
  if (defined($data[0])) {

    if ($data[0] eq "encoding:") {
      if (!defined($data[1]) || $data[1] ne "ASCII") {
        die "Error: ASCII encoding of NRRD file is required.";
      }
    }
    elsif ($data[0] eq "dimension:") {
      if (!defined($data[1])) {
        die "Error: Missing dimension.";
      }

      $dimension = $data[1];
    }
    elsif ($data[0] eq "sizes:") {
      shift @data;

      if (!defined($data[0]))
        { die "Error: Missing size information."; }

      @size = @data;
    }
  }
}

if (!defined($dimension)) { die "Error: Missing dimension."; };

if ($#size + 1 < $dimension) {
  my $axis = $#size+1;
  die "Error: Missing axis size information for axis $axis.";
}
if ($#size + 1 > $dimension) {
  my $axis = $#size+1;
  print "Warning:  Number of axis sizes exceeds volume dimension.\n";
  print "          Ignoring extra extra sizes.\n\n";
}

my $numv = 1;
for (my $d = 0; $d < $dimension; $d++)
  { $numv = $size[$d]*$numv; }


while (<INPUT>) {

  chomp;
  my $input_line = $_;
  $linenum++;

  my @data = split(" ",$input_line);

  push @scalar, @data;
}

close (INPUT);

if ($numv != $#scalar+1)
{  die "Error. Number of scalar values does not match array size."; };

open (OUTPUT, ">$output_filename") || 
    die ("$argv[1] cannot be opened for writing: $!\n");

foreach my $output_line (@header) {
  print (OUTPUT $output_line);
}

if (defined($flag_vector)) {
  my $k = 0;
  foreach my $s (@scalar) {
    print (OUTPUT $s);

    $k++;
    if (defined($size[0])) {
      if ($k%$size[0] == 0)
        { 
          print (OUTPUT "   "); 
          if (defined($size[1])) {
            if ($k%($size[0]*$size[1]) == 0)
              { print (OUTPUT "\n"); };

            if (defined($size[2])) {
              if ($k%($size[0]*$size[1]*$size[2]) == 0)
                { print (OUTPUT "\n"); };
            }
          }
        }
      else
        { print (OUTPUT " "); }
    }
  }

}
else {
  my $k = 0;
  foreach my $s (@scalar) {
    print (OUTPUT $s);

    $k++;
    if (defined($size[0])) {
      if ($k%$size[0] == 0)
        { 
          print (OUTPUT "\n"); 
          if (defined($size[1])) {
            if ($k%($size[0]*$size[1]) == 0)
              { print (OUTPUT "\n"); };
          }
        }
      else
        { print (OUTPUT " "); }
    }
  }
}

close(OUTPUT);



#---------------------------------------------------

sub parse_args {

  if (!defined($argv[0])) { &usage_error; };

  if ("$argv[0]" eq "-vector") {
    $flag_vector = 1; 
    shift @argv;
  }

  if (!defined($argv[0])) { &usage_error; };
  $input_filename = $argv[0];

  if (defined($argv[1]))
    { $output_filename = $argv[1]; }
  else
    { $output_filename = $argv[0]; }
}


#---------------------------------------------------
sub usage_error {
    print "Usage: [-vector] $0 {.nrrd ASCII file} [output .off file]\n";	
    exit;
};

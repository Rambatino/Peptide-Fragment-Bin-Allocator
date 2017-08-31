#!/usr/bin/perl
#
#
#
use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use vars qw($l $u $bin $startWith $contains $outPutFileName $empty $mA);
use POSIX;
use Math::Round 'nlowmult';

=comment
 All the variables from the command line:
 
        l => the lowest m/z value which is itself included in the search criteria (global variable = $lowerBound).
        u => the upper m/z value which is itself not included in the search criteria (global variable = $upperBound).
        b => the bin size (global variable = $binSize):
                
                                   b*(n) =< m/z < b*(n+1)
 
        mA => mass accuracy of the bins. The value must be positive and less than the
              bin size otherwise it will revert to zero and a fix bin output will be 
              created (global variable = $massAcc).
 
        s => one or more 'single-letter' amino acids, defining the starting sequence of included peptide
             chains (global variable = $startsWith).
        s => one or more 'single-letter' amino acids, defining part of the sequence of included peptide
             chains (global variable = $contains).
        o => string defining the name of the file (global variable = $outPutFileName)
 
=cut


Getopt::Long::Configure( qw(pass_through) );
GetOptions ('l=s' => \$l, 'u=s' => \$u, 'b=s' => \$bin, 'mA=s' =>\$mA, 's=s' => \$startWith, 'c=s' => \$contains, 'o=s' => \$outPutFileName);

#defining default values for lower bound, upper bound and the bin size if they have not been defined on
#the command line

my $lowerBound = ( defined $l ) ? $l : 1500;
my $upperBound = ( defined $u ) ? $u : 2000;
my $binSize = ( defined $bin ) ? $bin : 1;

#checking whether mass accuracy has been defined on the command line, and if it has whether it satisfies
#two prerequisites such that it is less than the bin size and also that it is greater than zero, otherwise
#it will default to zero and fixed bin sizes will be produced


my $massAcc = ( defined $mA ) ? $mA : $binSize;

#defining a global hash which is to be passed into the subroutines for modification

my %globalHash;

#
#
# Opening the file as defined as the first argument on the command line
# if it isn't the first argument, the program initiates a forced crash with a specific warning
# message
#
#

my $file = $ARGV[0];
    open(FILE,$file) or die "error: unable to open file $file\n";
    while ( <FILE> ) {
        my $temp = $_;
        my $dec_nums;
        my $aminoAcidString;

        #searching through each line and selecting the correct strings
        
        if($temp =~ /\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
            $dec_nums = $1;
            $aminoAcidString = $2;
        }
        
        #check to see whether the m/z value is a floating point number, as there is only one in the
        #file, if it isn't, $dec_nums is found by searching specifically for the floating point number
        
        #$dec_nums = temp =~ m{(\d+\.\d+)} if !isfloat $dec_nums;
        
        $dec_nums = "" unless defined $dec_nums; # make sure that string is not undefined
        if($dec_nums ne "" && $dec_nums >= $lowerBound && $dec_nums < $upperBound) { #initial check to determine whether the value falls within the upper and lower bound
            if ((defined $startWith) && substr($aminoAcidString, 0, length $startWith) eq $startWith) { #check to see if sequence starts with amino acids defined on command line
                addFloatTo($dec_nums,\%globalHash);
            } elsif ((defined $contains) && index($aminoAcidString, $contains) > 0) { #check to see if sequence contains amino acids as defined on command line
                addFloatTo($dec_nums,\%globalHash);
            } elsif(!defined $aminoAcidString){ #if there is no amino acid sequence found, then print message
                print "\nRecord found where no sequence of amino acid given\n";
            } else { #if neither sequence selection criteria have been defined, then all sequences within bounds are included
                addFloatToBin($dec_nums,\%globalHash);
        }
    }
}
writeHashToCSV("Bin","Number of Peptides",\%globalHash);

exit(1);

#######################################################################

######################## SUB ROUTINES #################################

#######################################################################



######

# Routine 1

=comment
 
 'addFloatToBin' subroutine takes the m/z floating point
 number as the first argument and a reference to the global hash as the
 second argument. It determines the bin(s) (i.e. key(s)) that the floating
 point number belongs to. The formula for the bin is:
 
 
 $binSize*floor(($float-$lowerBound)/$binSize) + $lowerBound
 + ($massAcc*floor(($float - $lowerBound)/$massAcc) + $lowerBound -($binSize*floor(($float-$lowerBound)/$binSize) + $lowerBound))
 - $count;
 
 
 The formula can be split into two parts:
 
        The first line determines the correct bin for the m/z value which is
 dependent upon the bin size and, if it is a bin that doesn't divide exactly
 into one, the number of bins that proceed it.
 
        The second two lines determine the error distribution around the m/z
 value and thus allow the value to be put into 'sliding bins' if a mass accuracy
 value is given. They calculate the size of the m/z value which is greater than
 the big bin value it falls into, and then loops through the while loop determining
 the smaller bins that the m/z value would fall into. Once again, the smaller bins
 are dependent upon the number of smaller bins that came before them and thus it has
 to calculate this as well (in the second floor() method).
 
 When no if no mass accuracy value is given, then it assumes that
 the mass accuracy is equal to the bin size and therefore the second two terms
 become = 0 and the while loop only completes its cycle once and thus forming a 
 fixed bin.

 For the PPM to daltons conversion, it assumes the bin sizes increase in proportion
 to the m/z value supplied and the sliding bins are given a generic '0.001' width. The
 conversion only occurs when the mass accuracy value exceeds that of the bin size.
 
 It uses this formula for three reasons:
 
 a) bin sizes that aren't divisible by one, give un-even bin estimates
 when finding the lower bound of the floating point number when using
 nlowmult(). eg. a 0.3 bin size may yield a bin estimate as 1511.2 =<
 m/z < 1511.5 (based upon the proceeding number of bins) which wouldn't
 be reproducible by setting the lower bound estimate to be 0.3 using
 nlowmult() as is the same in the second part of the equation.
 
 b) defining the bins/keys on the fly is much faster than creating all the
 bins/keys seperately, searching through each of the keys of the hash
 and then incrementing their values. The formula dictates the keys and
 therefore it is O(N) as opposed to O(N^N).
 
 c) it means that only one sub routine is required for sliding and fixed bins
 as well as tandemly allow switching between PPM and Da.
 
 
=cut


sub addFloatToBin {
    my ($float,$originalHash) = @_;
    my $count = 0;
    if ($massAcc > $binSize) {
        $binSize = PPMtoDaltons($massAcc,$float);
        $massAcc = 0.001;
    }
    while ($count < $binSize) {
        my $key = $binSize*floor(($float-$lowerBound)/$binSize) + $lowerBound + ($massAcc*floor(($float - $lowerBound)/$massAcc) + $lowerBound-($binSize*floor(($float-$lowerBound)/$binSize) + $lowerBound)) - $count;
        if(exists $originalHash->{$key}){
            ++$originalHash->{$key};
        } else {
            $originalHash->{$key} = 1;
        }
        $count += $massAcc;
    }
}


# Routine 4

=comment
 
 'writeHashToCSV' subroutine takes three arguments:
        The name of the first column;
        The name of the second column;
        The hash in which the keys fill the first column with their
        associated values in the second column.
 It then produces a .csv file. The name is alterable on the command line.
 
=cut


sub writeHashToCSV {
    my ($firstColumnHeader, $secondColumnValues, $hash) = @_;
    my $filename = (defined $outPutFileName)?sprintf('%s.csv',$outPutFileName):'histogramData.csv';
    open(my $fh, ">", $filename) or die "Could not open output file '$filename' $!";
    print $fh "$firstColumnHeader,$secondColumnValues\n";
    while(my($k, $v) = each %{$hash}) {
            print $fh "$k,$v\n";
    }
    close $fh;
}

sub PPMtoDaltons {
    #takes two arguments, the first is the mass accuracy in parts per million
    #and the second is the m/z that the error is about
    return ($_[0] * $_[1])/(10**6);
}



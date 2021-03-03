#!/usr/bin/perl
# LOPT: program for Level OPTimization
# Author: A. Kramida
# Version 2.3
# Changes compared to v. 2.1: Minor bug fixes.
#   1) In ASCII output for lines, sometimes some columns were concatenated.
#      This was fixed in v. 2.1.
#   2) If in the input lines file the level labels were short, e.g.,
#      numerical idexes were used, the program falsely reported duplicate
#      transitions. This was fixed in v. 2.2 by inserting a double-comma separator in
#      the search key of the line list sorted by level names.
#   3) In the include file vacair.pl, the tolerance value (variable $tol)
#      was increased in v. 2.2 from 1e-19 to 1e-12 due to insufficient numerical precision
#      of floating-point arithmetics in Perl, which sometimes caused the
#      program to hang up while converting air wavelengths to vacuum.
#   4) An error in the Round() function has been corrected in v. 2.3 by
#      M. Ya Vokhmentsev. This error was causing an incorrect rounding of
#      negative numbers.
#   5) Errors in the code that caused wrong results when working with
#      variable substitution with multiple fixed levels,and/or division
#      into independent groups have been fixed in v. 2.3 by A. Kramida.
#   6) Error in handling cases involving transitions with the ground level
#      on the right side (with negative observed energy) fixed in v. 2.3
#      by A. Kramida
# History:
# Version 2.1 (first Perl version) released in 2006 and published in CPC.
# Version 2.2 of 2013 (unpublished)
# Version 2.3, August 2017 (unpublished)
#
use strict;
require 'InvertMatrix.pl';
require 'vacair.pl';

use FindBin qw($RealBin);
use lib $RealBin;

my $LevRoundThreshold = 25.0;
my $TrRoundThreshold  = 25.0;
my $MinAirWn = 5000.0; # cm-1
my $MaxAirWn = 50000.0; # cm-1
my $ParName = 'lopt.par'; # Default par. file name

     # Flag masks
my $Fair          =    1;
my $FBlend        =    2;
my $FQuest        =    4;
my $FMask         =    8;
my $FPredict      =   16;
my $FExcludedTo   =  128;
my $FExcludedFrom =  256;
my $FVirtual      =  512;
my $FExcluded     = 1024;
my $FMP           =   24; # Masked OR Predicted

my $too_small     = 1e-40;
my $eps = &get_eps;

# Wavelength->energy conversion factors; NIST Fund. Const. (CODATA 2006)
# Negative values indicate that the conversion involves reciprocation
my %WLfactors =
 ('A'   => 1.0,
  'nm'  => 0.1,
  'mcm' => 1e-4,
  'cm-1'=>-1e8,
  'eV'  =>-1.239841875e4,
  'MHz' =>-2.997924580e12,
  'GHz' =>-29.97924580e8);

#      1 A = 1e8/cm-1 = 1e8/2.997924580e4 MHz (exact)
#      1 cm-1 = 1.239841875e-4 eV
#      conversion factors are taken from
#      2006 CODATA recomm}ed values,
#      http://physics.nist.gov/cuu/Constants/index.html

        # Fixing levels is done by assigning uncertainty 1/sqrt($FixLevWeight)
        # to all 'fixed' levels. Its value is multiplied by $FixWeightFactor.
        # Prior to v. 2.3, its value was set to $FixWeightFactor times the
        # smallest uncertainty of observed transitions and fixed levels
#my $FixedLevUncert;
#my $FixedLevWeight;
my $FixWeightFactor = 5e-4; # Weight factor used to effectively fix the "fixed" levels.

# Global variables
my (
  @E,                # array of level energies or uncertainties
  $A,                # array of matrix coefficients
  $B,                # inverse matrix
  $Fa,$fb,$Fax,      # Matrix Fa and vector fb containing estimates of floating-point errors of A and E; Fax = Fa*x
  %LinesByW,         # hash keyed by observed wavelength/wavenumber
  %LinesByLevels,    # hash keyed by lower/upper level names
  %Levels);          # hash keyed by level label
my (@Quest1,@Quest2, @Sigma, @FixLevs,@Lines,@Levels,@GroupLines,@CGroups,
  @UncBases,@SLLines,$SkipLambdas,$TuneLevels,$PrintUncCorr,$Uncert,$InvChk,
  $RoundInput,$WriteVirt,$TabDelimitedOutput,$NLines,$LineNo,$NA,$IntFirstCol,
  $IntLastCol,$WLFirstCol,$WLLastCol,$widthInt,$WuncFirstCol,$WuncLastCol,
  $WeightFirstCol,$WeightLastCol,$LLFirstCol,$LLLastCol,$LNameWidth,
  $ULFirstCol,$FlgFirstCol,$WunitFirstCol,$Groups,$Ground,$PLg,$MinUnc,
  $DivideGroups,$numFixLevs,$LastWunit,$LastWLUnc,$WLfactor,$CUnc,
  $max_num_error_frac,$ksi,$r,$cond,$rel_cond,$iter_crit,$use_substitution,
  $normwise_error_bound,$statistics_size,$trial_size,$use_Cholesky,
  $TransFilename,$FixLevFilename,$LevOutFilename,$TransOutFilename,@VirtLines);
my $MaxCGNo = 0;
my %Widths = ();

# General-purpose functions
############################################################################
sub get_eps {  #11/18/2009 10:10AM
############################################################################
  # Determine the machine floating-point precision, or machine epsilon
  my $x = 2;
  my $e = 1;
  while ( $x - 1 > 0 ) {
    $e /= 2;
    $x = 1 + $e;
  }
  return $e*2;
} ##get_eps

sub AirWavenumber($) {
  my $E = abs(shift);
  return (($E >= $MinAirWn) && ($E <= $MaxAirWn));
}

############################################################################
sub InitLevel() {   #10/7/2006 5:16PM A.Kramida
############################################################################
  # Initialization of a data structure for the Energy Level entity
  my $L = {};
  $L->{'Name'} = '';
  $L->{'To'}   = [];
  $L->{'From'} = [];
  $L->{'C'}    = [];
  $L->{'Subst'}    = undef;
  foreach (qw{Fixed GNo E D D1 D2 NT NT1 NT2 NB NL NQ NM ND NV Er Shift Precision d d1 Dr D1r D2r dr} ) {
    $L->{$_} = 0;
  }
  return $L;
} ##InitLevel()

############################################################################
sub InitTrans() {  #10/7/2006 5:22PM A.Kramida
############################################################################
  # Initialization of a data structure for the Transition entity
  my $T = {};
  $T->{'FromLevel'} = undef;
  $T->{'ToLevel'}   = undef;
  $T->{'Wunit'}    = 'A';
  $T->{'Weight'}    = 1;
  $T->{'C'}         = [];

  foreach (qw{Wobs Eobs WnUnc WnUncCalc Flags InGroup NumAssign CGroupNo ObsPrecision CalcPrecision}) {
    $T->{$_} = 0;
  }
  foreach (qw{_Lobs _Eobs _LuncObs _wnObs _WnUncObs _LairCalc _LvacCalc
    _LuncCalc _WnCalc _WnUncCalc _dW _dE})
  {
    $T->{$_} = '_';
  }
  foreach (qw{Iobs _flags}) {
    $T->{$_} = '';
  }
  return $T;
} ##InitTrans

############################################################################
sub trans_Ecalc($) {
############################################################################
  # Return the Ritz wavenumber of a transition (unrounded).
  # The {'Shift'} field holds the adjustment of energy
  # that may be necessary because of single-line levels.
  my $T = shift;
  return $T->{'FromLevel'}->{'E'} + $T->{'FromLevel'}->{'Shift'} - $T->{'ToLevel'}->{'E'} - $T->{'ToLevel'}->{'Shift'};
}

############################################################################
sub trans_Ecalc_from_rounded_levels($) {
############################################################################
  # Return the (unrounded) Ritz wavenumber of a transition, as calculated
  # from rounded level values
  my $T = shift;
  return $T->{'FromLevel'}->{'Er'} - $T->{'ToLevel'}->{'Er'};
}

############################################################################
sub trans_Ecalc_rounded($) {
############################################################################
  # Return the rounded Ritz wavenumber of a transition, as calculated
  # from rounded level values
  my $T = shift;
  my $E = $T->{'FromLevel'}->{'Er'} - $T->{'ToLevel'}->{'Er'};
  my $Eunc = ($Uncert ? $T->{'WnUncCalc'} : 0.01);
  my ($Xr,$u) = round($E,$Eunc,$TrRoundThreshold,1);
  $u = '_' unless $Uncert;
  return ($Xr,$u);
}

############################################################################
sub trans_air($) {   #11/11/2009 9:50AM
############################################################################
  # Return 1 if wavelength of a transition is in air, 0 otherwize
  my $T = shift;
  return (($T->{'Flags'} & $Fair) || ($T->{'Eobs'}+0 && AirWavenumber($T->{'Eobs'})) || AirWavenumber(trans_Ecalc($T)));
} ##trans_air($)

############################################################################
sub trans_WLvac_calc($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the vacuum Ritz wavelength of a transition (unrounded)
  my $T = shift;
  my $Wn = trans_Ecalc($T);
  return (abs($Wn) < $too_small ? 0 : 1e8/$Wn);
}

############################################################################
sub trans_WLvac_calc_rounded($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the rounded vacuum Ritz wavelength of a transition
  # as calculated from the rounded level values.
  my $T = shift;
  my ($Wn,$WnUnc) = trans_Ecalc_from_rounded_levels($T);
  return ('_','_') if (abs($Wn) < $too_small);
  $WnUnc = $T->{'WnUncCalc'};
  my $lambda = 1e8/$Wn;
    # uncertainty of the predicted wavelength
  my $WLunc = ( $Uncert ? abs(($WnUnc/$Wn)*$lambda) : 0.001 );
  my ($Xr,$u) = round($lambda,$WLunc,$TrRoundThreshold,1);
  $u = '_' unless $Uncert;
  return ($Xr,$u);
}

############################################################################
sub trans_WLair_calc_rounded($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the rounded air Ritz wavelength of a transition
  # as calculated from the rounded level values.
  my $T = shift;
  my ($Wn,$WnUnc) = trans_Ecalc_from_rounded_levels($T);
  return ('_','_') if (abs($Wn) < $too_small);
  $WnUnc = $T->{'WnUncCalc'};
  my $lambda = Lair(1e8/abs($Wn));
  $lambda = -$lambda if ($Wn < 0);
    # uncertainty of the predicted wavelength
  my $WLunc = ( $Uncert ? abs(($WnUnc/$Wn)*$lambda) : 0.001 );
  my ($Xr,$u) = round($lambda,$WLunc,$TrRoundThreshold,1);
  $u = '_' unless $Uncert;
  return ($Xr,$u);
}

############################################################################
sub trans_get_WL_from_Wn($$) {   #11/11/2009 9:52AM
############################################################################
  # Return an air or vacuum wavelength corresponding to the given wavenumber
  # (second parameter), based on the properties of a given transition
  # (first parameter).
  my ($T,$Wn) = @_;
  return 0 if (abs($Wn) < $too_small);
  my $WL;
  if ( trans_air($T) ) {
    $WL = Lair(1e8/abs($Wn));
    $WL = -$WL if ($Wn < 0);
  } else {
    $WL = 1e8/$Wn;
  }
  return $WL;
} ##trans_get_WL_from_Wn($$)

############################################################################
sub trans_WLobs($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the wavelength corresponding to the transition's observed
  # wavenumber. The wavelength will be in air or in vacuum, depending on
  # the transition.
  my $T = shift;
  return trans_get_WL_from_Wn($T, $T->{'Eobs'});
} ##trans_WLobs($)

############################################################################
sub trans_WLobs_rounded($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the rounded wavelength corresponding to the transition's observed
  # wavenumber. The wavelength will be in air or in vacuum, depending on
  # the transition.
  my $T = shift;
  my $WL = trans_get_WL_from_Wn($T, $T->{'Eobs'});
  return ('_','_') unless ( $WL && !($T->{'Flags'} & $FVirtual) );
  my $WLunc = $T->{'WnUnc'}/$T->{'Eobs'}*$WL;
  my ($Xr,$u) = round($WL,$WLunc,$TrRoundThreshold,1);
  return ($Xr,$u);
} ##trans_WLobs_rounded($)

############################################################################
sub trans_WnObs($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the observed wavenumber of the transition.
  return shift->{'Eobs'};
} ##trans_WnObs($)

############################################################################
sub trans_WnObs_rounded($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the rounded observed wavenumber of the transition.
  my $T = shift;
  return ('_','_') if ($T->{'Flags'} & $FVirtual);
  return ('_','_') if ( ($T->{'Flags'} & $FPredict) && !$T->{'Eobs'} );
  my $Wn = $T->{'Eobs'};
  my $WnUnc = $T->{'WnUnc'};
  my ($Xr,$u) = round($Wn,$WnUnc,$TrRoundThreshold,1);
  return ($Xr,$u);
} ##trans_WnObs_rounded($)

############################################################################
sub trans_WLcalc($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the (unrounded) Ritz wavelength of a transition (air or vacuum).
  my $T = shift;
  return trans_get_WL_from_Wn($T, trans_Ecalc($T));
} ##trans_WLcalc($)

############################################################################
sub trans_WLcalc_from_rounded_levels($) {   #10/7/2006 6:09PM A.Kramida
############################################################################
  # Return the (unrounded) Ritz wavelength of a transition (air or vacuum),
  # calculated from rounded level values.
  my $T = shift;
  return trans_get_WL_from_Wn($T, trans_Ecalc_from_rounded_levels($T));
} ##trans_WLcalc($)

############################################################################
sub trans_WLcalc_rounded($) {
############################################################################
  # Return the rounded Ritz wavelength of a transition (air or vacuum),
  # calculated from rounded level values.
  my $T = shift;
  if ( trans_air($T) ) {
    return trans_WLair_calc_rounded($T);
  } else {
    return trans_WLvac_calc_rounded($T);
  }
}

sub trans_DoublyAssigned($) {
  # Return 1 if a spectral line corresponding to the transition is
  # multiply assigned (assigned to several different transitions).
  my $T = shift;
  return $T->{'NumAssign'} > 1;
}

sub trans_SingleLineLevelFlag($) {
  # Return star ('*') if one of the levels involved in the transition is
  # determined solely by this transition.
  my $T = shift;
  my $SingleLineLevelFlag =
  (!($T->{'Flags'} & ($FMP | $FVirtual)) &&
    (($T->{'FromLevel'}->{'NT2'} + $T->{'FromLevel'}->{'NV'} == 1)
     || ($T->{'ToLevel'}->{'NT2'} + $T->{'ToLevel'}->{'NV'} == 1))
  ) ? '*' : ' ';
  return $SingleLineLevelFlag;
}

############################################################################
sub trans_ObsQuantityValue($) {  #11/11/2009 12:55PM
############################################################################
  # Return the value of the observed quantity as read from the transitions
  # input file.
  my $T = shift;
  return $T->{'Wobs'};
} ##trans_ObsQuantityValue($)

############################################################################
sub trans_ObsQuantityUnc($) {  #11/11/2009 12:55PM
############################################################################
  # Return the uncertainty of the observed quantity as read from the
  # transitions input file.
  my $T = shift;
  if ( $T->{'Flags'} & ($FMP | $FVirtual) ) {
    return '_';
  }
  my $WnUnc = $T->{'WnUnc'};
  my $WLfactor = $WLfactors{$T->{'Wunit'}};
  my $u;
  if ($WLfactor < 0) {
      # Energy units
    $u = -1e-8*$WLfactor*$WnUnc;
  } else {
      # Angstroms  or micrometers
    my $Wn = $T->{'Eobs'};
    my $W = $WLfactor*trans_WLobs($T);
    $u = $WnUnc/$Wn*$W;
  }
  if ( $u < $too_small) {
    $u = '_';
  }
  return $u;
} ##trans_ObsQuantityUnc($)

############################################################################
sub trans_ObsQuantityRounded($) {   #11/11/2009 12:53PM
############################################################################
  # Return the rounded value of the observed quantity.
  my $T = shift;
  my $W = trans_ObsQuantityValue($T);
  if ( $T->{'Flags'} & ($FMP | $FVirtual) ) {
    return ('_','_');
  }
  my $Wunc = trans_ObsQuantityUnc($T);
  my ($Xr,$u) = round($W,$Wunc,$TrRoundThreshold,1);
  return ($Xr,$u);
} ##trans_ObsQuantityRounded($)

############################################################################
sub trans_CalcValueOfObsQuantity($) {  #11/10/2009 12:56PM
############################################################################
  # Return the Ritz value of the observed quantity W in the same units as W.
  my $T = shift;
  my $WLfactor = $WLfactors{$T->{'Wunit'}};
  my $Wc;
  if ($WLfactor < 0) {
      # Energy units
    $Wc = -1e-8*$WLfactor*trans_Ecalc($T);
  } else {
      # Angstroms  or micrometers
    $Wc = $WLfactor*trans_WLcalc($T);
  }
  return $Wc;
} ##trans_CalcValueOfObsQuantity($)

############################################################################
sub trans_CalcValueOfObsQuantityFromRoundedLevels($) {  #11/10/2009 12:56PM
############################################################################
  # Return the Ritz value of the observed quantity W in the same units as W,
  # calculated from rounded level values.
  my $T = shift;
  my $WLfactor = $WLfactors{$T->{'Wunit'}};
  my $Wc;
  if ($WLfactor < 0) {
      # Energy units
    $Wc = -1e-8*$WLfactor*trans_Ecalc_from_rounded_levels($T);
  } else {
      # Angstroms  or micrometers
    $Wc = $WLfactor*trans_WLcalc_from_rounded_levels($T);
  }
  return $Wc;
} ##trans_CalcValueOfObsQuantity($)

############################################################################
sub trans_CalcValueOfObsQuantityUnc($) {  #11/11/2009 12:55PM
############################################################################
  # Return the uncertainty of the Ritz value of the observed quantity W
  # in the same units as W.
  my $T = shift;
  my $WnUnc = $T->{'WnUncCalc'};
  my $WLfactor = $WLfactors{$T->{'Wunit'}};
  my $u;
  if ($WLfactor < 0) {
      # Energy units
    $u = -1e-8*$WLfactor*$WnUnc;
  } else {
      # Angstroms  or micrometers
    my $Wn = trans_Ecalc($T);
    my $W = $WLfactor*trans_WLcalc($T);
    $u = $WnUnc/$Wn*$W;
  }
  return $u;
} ##trans_CalcValueOfObsQuantityUnc($)

############################################################################
sub trans_CalcValueOfObsQuantityRounded($) {  #11/10/2009 12:56PM
############################################################################
  my $T = shift;
  my $Wc = trans_CalcValueOfObsQuantityFromRoundedLevels($T);
  my $WcUnc = trans_CalcValueOfObsQuantityUnc($T);
  my ($Xr,$u,$lsd) = round($Wc,$WcUnc,$TrRoundThreshold,1);
  return ($Xr,$u,$lsd);
} ##trans_CalcValueOfObsQuantityRounded($)

############################################################################
sub straight_sum($) {  #12/14/2009 12:32PM
############################################################################
  # Calculate the straight sequential sum of a given set of values.
  # This summation method has bad numerical properties and is not used in
  # the program. The function is provided for testing purposes only.
  my $v = shift;
  my $N = $#{$v} + 1;
  return 0 unless $N;
  my $S = 0;
  for ( my $j = 0; $j < $N; $j++ ) {
    $S += $v->[$j];
  }
  return $S;
  # Mean numerical error of straight summation of N terms has been found to
  # obey e = 0.148*sqrt(N)*sqrt(sum(x_i*x_i))*eps
} ##straight_sum($)

############################################################################
sub Kohan_sum($) {   #12/14/2009 11:54AM
############################################################################
  # Return the sum of a given set of values, calulated using the Kohan method.
  my $v = shift;
  my $N = $#{$v} + 1;
  return 0 unless $N;
  my $S = $v->[0];
  my $C = 0;
  for ( my $j = 1; $j < $N; $j++ ) {
    my $Y = $v->[$j] - $C;
    my $T = $S + $Y;
    $C = ($T - $S) - $Y;
    $S = $T;
  }
  # Mean numerical error of Kohan sum of N terms has been found to
  # obey e = 0.27*sqrt(sum(x_i*x_i))*eps (see below).
  return $S;
} ##Kohan_sum($)

############################################################################
sub Kohan_sum_error($) {  #12/14/2009 11:59AM
############################################################################
  # Return an estimate of the numerical error in Kohan summation of a given
  # set of floating-point numbers.
  my $v = shift;
  my $N = $#{$v} + 1;
  my $S = 0;
  return $S unless $N;
  for ( my $j = 0; $j < $N; $j++ ) {
    my $x = $v->[$j];
    $S += $x*$x;
  }
  # The factor 0.27 results from numericval experiments comparing randomly
  # generated sums calculated in real*4 and real*8 formats.
  # It has been found to hold at least within interval (2,100000) for
  # the number of terms in the sum.
  $S = 0.27*$eps*sqrt($S);
  return $S;
} ##Kohan_sum_error

sub Trim($) {
  # Trim leading and trailing blanks from the given string.
  # If the string starts with a period, add a zero before this period
  # (to ensure proper treatment of real floating-point numbers
  # given without the leading zero).

  my $S = shift;
  $S =~ s/^\s+|\s+$//g;
  $S =~ s/^[.]/0./;
  return $S;
}

sub Error($) {
  # Report an error and terminate the process
  print "\n",shift,"\n";
  exit(1);
}

sub MergeGroups($$$;$) {
  # Merge two groups of levels, $G1 and $G2, from the list of groups $G,
  # into one group; replace the two entries in the group list with the
  # resulting merged group.
  my ($G1,$G2,$Groups,$FixHash) = @_;

  # Select the largest group as $G
  my ($G, $G3) = ($#{$G1} > $#{$G2}) ? ($G1, $G2) : ($G2, $G1);

  my $k = $#{$G};

  # Insert elements (levels) of group $G3 into group $G
  push(@{$G},@{$G3});

  if ( $FixHash ) {
    foreach my $PL (@{$G3}) {
      #push(@{$G},$PL);

      $k++;
      $Levels{$PL->{'Name'}} = [$G,$k];
    }
  }

  # Delete group $G3 from the list of groups
  for ( my $i = 0; $i <= $#{$Groups}; $i++ ) {
    if ( $Groups->[$i] eq $G3 ) {
      splice(@{$Groups}, $i, 1);
      last;
    }
  }
  return $G;
}

sub SameLevels($$) {
  # return 1 if the two levels' labels are the same
  my ($PLx,$PLy) = @_;
  return ($PLy->{'Name'} eq $PLx->{'Name'});
}

sub InsertLevs($) {
  # Insert levels of a given transition into the list of levels
  my $PT = shift;

  my ($G1,$G2,$PL1,$PL2,$PT1) = (undef,undef,undef,undef,undef);
  my $PL = $PT->{'FromLevel'};

  # Try to find lower and upper levels in already existing lists
  my ($k1,$k2);
  ($G1,$k1) = defined($Levels{$PL->{'Name'}}) ? @{$Levels{$PL->{'Name'}}} : (undef,-1);
  if ( defined($G1) ) {
    $PL1 = $G1->[$k1];
    if ($PL1 ne $PL) {
      $PT->{'FromLevel'} = $PL1;
    }
  }

  $PL = $PT->{'ToLevel'};
  ($G2,$k2) = defined($Levels{$PL->{'Name'}}) ? @{$Levels{$PL->{'Name'}}} : (undef,-1);
  if ( defined($G2) ) {
    $PL2 = $G2->[$k2];
    if ($PL2 ne $PL) {
      $PT->{'ToLevel'} = $PL2;
    }
  }

  # Keep line counters for both levels
  unless ($PT->{'Flags'} & $FMP) {
    if ($PT->{'Flags'} & $FQuest) {
      $PT->{'FromLevel'}->{'NQ'}++;
      $PT->{'ToLevel'}->{'NQ'}++;
    } elsif ($PT->{'Flags'} & $FBlend) {
      $PT->{'FromLevel'}->{'NB'}++;
      $PT->{'ToLevel'}->{'NB'}++;
    }
  }
  if ($G1) {
    if ($G2) {
      if (!($PT->{'Flags'} & $FVirtual)) {
        my $i_prev = $LinesByLevels{$PT->{'ToLevel'}->{'Name'} . ',,' . $PT->{'FromLevel'}->{'Name'}};
        if ($i_prev) {
          print "Warning: duplicate transitions at lines no. $LineNo and $i_prev\n";
        }
        $i_prev = $LinesByLevels{$PT->{'FromLevel'}->{'Name'} . ',,' . $PT->{'ToLevel'}->{'Name'}};
        if ($i_prev) {
          print "Warning: duplicate reverse transitions at lines no. $LineNo and $i_prev\n";
        }
      }
      if (($G1 ne $G2) and !($PT->{'Flags'} & $FMP)) {
        $G1 = MergeGroups($G1,$G2,$Groups,1);
      }
    } else {
      $PL2 = $PT->{'ToLevel'};
      push(@{$G1}, $PL2);
      $Levels{$PL2->{'Name'}} = [$G1,$#{$G1}];
      push(@Levels,$PL2);
    }
  } else {
    if ($G2) {
      $PL1 = $PT->{'FromLevel'};
      push(@{$G2},$PL1);
      $Levels{$PL1->{'Name'}} = [$G2,$#{$G2}];
      push(@Levels,$PL1);
    } else {
      $PL1 = $PT->{'FromLevel'};
      $PL2 = $PT->{'ToLevel'};
      $G1 = [$PL2];
      $Levels{$PL2->{'Name'}} = [$G1,$#{$G1}];
      push(@{$Groups},$G1);
      push(@Levels,$PL2);
      push(@Levels,$PL1);
      if (!($PT->{'Flags'} & $FMP)) {
        push(@{$G1},$PL1);
        $Levels{$PL1->{'Name'}} = [$G1,$#{$G1}];
      } else {
        $G2 = [$PL1];
        push(@{$Groups},$G2);
        $Levels{$PL1->{'Name'}} = [$G2,$#{$G2}];
      }
    }
  }
}

sub sqr($) {
  # Fail-safe power of 2 function
  my $x = abs(shift);
  return 1.79769312951459e+308 if ($x > 1.340780791e+154);
  return $x*$x;
}

############################################################################
sub ValidNumber($) {   #11/04/2006 10:47AM A.Kramida
############################################################################
  # Return 1 if the given string is a valid number, 0 otherwise.
  my $s = shift;
  return 0 if ( $s !~ /^[0-9.e+-]+$/i );
  my $n = $s + 0;
  return 0 if ( !$n && ($s !~ /^0*(\.0+){0,1}$/) );
  return 1;
} ##ValidNumber($)

# Hash defining the valid parameter options
my %ParamDescription =
('TRANSITIONS INPUT FILE NAME'              => ['$TransFilename','filename','<','trans-inp.txt'],
 'FIXED LEVELS INPUT FILE NAME'             => ['$FixLevFilename','filename','<','fixlev-inp.txt'],
 'LEVELS OUTPUT FILE NAME'                  => ['$LevOutFilename','filename','>','lev-out.txt'],
 'TRANSITIONS OUTPUT FILE NAME'             => ['$TransOutFilename','filename','>','trans-out.txt'],
 'OMIT CALC. WAVELENGTHS'                   => ['$SkipLambdas','bool','YN',0],
 'TUNE SINGLE-LINE LEVELS'                  => ['$TuneLevels','bool','YN',1],
 'PRINT LISTING OF CORRELATED-LINES UNCERT' => ['$PrintUncCorr','bool','YN',0],
 'CALCULATE PREDICTED LINE UNCERTAINTIES'   => ['$Uncert','bool','YN',1],
 'USE CHOLESKY'                             => ['$use_Cholesky','bool','YN',1],
 'USE VARIABLE SUBSTITUTION'                => ['$use_substitution','bool','YN',1],
 'NUMBER OF TRIALS'                         => ['$statistics_size','int','>-1',0],
 'PERFORM MATRIX-INVERSION CHECK'           => ['$InvChk','bool','YN',0],
 'WRITE VIRTUAL LINES'                      => ['$WriteVirt','bool','YN',0],
 'DIVIDE LEVELS INTO INDEPENDENT GROUPS'    => ['$DivideGroups','bool','YN',0],
 'MIN. WAVENUMBER FOR AIR WAVELENGTH'       => ['$MinAirWn','float',undef,5000],
 'MAX. WAVENUMBER FOR AIR WAVELENGTH'       => ['$MaxAirWn','float',undef,50000],
 'ROUND-OFF THRESHOLD FOR OUTPUT LEVELS'    => ['$LevRoundThreshold','float',[5,99.99999],25],
 'ROUND-OFF THRESHOLD FOR OUTPUT LINES'     => ['$TrRoundThreshold','float',[5,99.99999],25],
 '1ST COLUMN OF LINE INTENSITY'             => ['$IntFirstCol','int','>0',1000],
 'LAST COLUMN OF LINE INTENSITY'            => ['$IntLastCol','int','>0',1000],
 '1ST COLUMN OF WAVELENGTH VALUE'           => ['$WLFirstCol','int','>0',undef],
 'LAST COLUMN OF WAVELENGTH VALUE'          => ['$WLLastCol','int','>0',undef],
 '1ST COLUMN OF WAVELENGTH UNITS'           => ['$WunitFirstCol','int','>0',1000],
 '1ST COLUMN OF WAVELENGTH UNCERTAINTY'     => ['$WuncFirstCol','int','>0',undef],
 'LAST COLUMN OF WAVELENGTH UNCERTAINTY'    => ['$WuncLastCol','int','>0',undef],
 '1ST COLUMN OF LINE WEIGHT'                => ['$WeightFirstCol','int','>0',1000],
 'LAST COLUMN OF LINE WEIGHT'               => ['$WeightLastCol','int','>0',1000],
 '1ST COLUMN OF LOWER LEVEL LABEL'          => ['$LLFirstCol','int','>0',undef],
 'LAST COLUMN OF LOWER LEVEL LABEL'         => ['$LLLastCol','int','>0',undef],
 '1ST COLUMN OF UPPER LEVEL LABEL'          => ['$ULFirstCol','int','>0',undef],
 '1ST COLUMN FOR LINE FLAGS'                => ['$FlgFirstCol','int','>0',1000],
 'TAB-DELIMITED OUTPUT'                     => ['$TabDelimitedOutput','bool','YN',1],
 'ROUND INPUT'                              => ['$RoundInput','bool','YN',1]
);

############################################################################
sub AssignParam($$) {   #10/16/2009 11:27AM
############################################################################
  # Assign a given value (second parameter) to a global variable
  # defined by a given parameter key (first parameter).
  my ($key,$value) = @_;
  my ($var_name,$type,$valid_range,$default_value) = @{$ParamDescription{$key}};

  $value =~ s/^\s*([^;]+);.*$/$1/;   # delete everything after semicolon
  $value =~ s/^\s+|\s+$//g;          # Trim leading and trailing spaces

  $ParamDescription{$key}->[4] = 1;  # Mark the parameter as assigned

  $key = lc($key);
  my $use_default = 0;
  if ( $value eq '' ) {
    if ( defined($default_value) ) {
      $value = $default_value;
      $use_default = 1;
    } else {
      Error("Value for parameter \"$key\" is missing in the parameter file.")
    }
  }
  my $given = $use_default ? "assumed as deafult" : "specified in the parameter file";
  my $val = uc($value);
  if ( $type eq 'filename' ) {
    my $file_mode = $valid_range;
    unless ( ($file_mode eq '>') || (-f $value) ) {
      Error("File $value for \"$key\" $given does not exist.");
    }
    my $dealt_with = ($file_mode eq '>' ? 'created' : 'opened');
    open(TMP, "$file_mode$value") or
      Error("File $value for \"$key\" $given could not be $dealt_with.");
    close TMP;
    $value = "qw($value)";
  } elsif ($type eq 'bool') {
    $val = 'Y' if $value eq '1';
    $val = 'N' if $value eq '0';
    unless ( $val =~ /^[YN]$/ ) {
      Error("Invalid value \"$value\" is specified for [arameter \"$key\". Valid values are \"Y\" or \"N\".");
    }
    $value = ($val eq 'Y') ? 1 : 0;
  } elsif ($type eq 'int') {
    $val =~ s/^0+([^0].*)$/$1/;              # drop leading zeros
    $val =~ s/^0+$/0/;                       # drop leading zeros
    $val = "0$val" if $val =~ /^\./;         # add a zero at the start if the value starts with a decimal point
    $val =~ s/(\.\d*[1-9])0*$/$1/;           # drop trailing zeros
    $val =~ s/\.0*$//;                       # drop trailing zeros and a decimal point
    $val = '0' if ($val =~ /^-0*$/);         # drop the minus sign before zero values
    if ( $val ne int($val) ) {
      Error("Invalid format of value \"$value\" for parameter \"$key\". It must be an integer number.");
    }
    $value = $val;
  } elsif ($type eq 'float') {
    $val =~ s/^0+([^0].*)$/$1/;              # drop leading zeros
    $val =~ s/^0+$/0/;                       # drop leading zeros
    $val = "0$val" if $val =~ /^\./;         # add a zero at the start if the value starts with a decimal point
    $val =~ s/(\.\d*[1-9])0*$/$1/;           # drop trailing zeros
    $val =~ s/\.0*$//;                       # drop trailing zeros and a decimal point
    if ( $val ne ($val + 0) ) {
      Error("Invalid format of value \"$value\" for parameter \"$key\". It must be a floating-point number.");
    }
    $value = $val;
  }
  if ((($type eq 'int') || ($type eq 'float')) && defined($valid_range) ) {
    my ($range_type, $min, $max);
    if ( $valid_range =~ /^[<>\[\]]([^\[\]]+)\]*$/ ) {
      ($range_type,$min) = ($1,$2);
      if ( (($range_type eq '>') && ($value <= $min)) || (($range_type eq '<') && ($value >= $min)) ) {
        Error("Value \"$value\" specified for parameter \"$key\" is outside the valid range. It must be $valid_range.");
      } elsif ($range_type eq '[') {
        $min =~ /^([^,]+),([^,+])/;
        ($min,$max) = ($1,$2);
        if ( ($value < $min) || ($value > $max) ) {
          Error("Value \"$value\" specified for parameter \"$key\" is outside the valid range. It must be in the interval $valid_range.");
        }
      }
    }
  }
  if ($var_name) {
    eval("$var_name=$value");
  } else {
    $ParamDescription{uc($key)}->[0] = $value;
  }
} ##AssignParam($$)

sub ReadParams($) {
  # Read parameter values from the input parameter file
  # and assign these values to the corresponding global variables.
  my $par = shift;
  $ParName = $par if $par;

  open(PARFILE, $ParName) or Error("Cannot open parameters file $ParName");

  my %Params = ();
  my $s;
  while (($s = <PARFILE>) && defined($s)) {
    my $key = undef;
    $s =~ s/^([^;]+);(.*)first(.*)$/$1;${2}1st$3/i;
    $s =~ s/^([^;]+);(.*)level name(.*)$/$1;${2}level label$3/i;
    $s =~ s/WAVELENGTH in/WAVELENGTH VALUE IN/i;
    foreach my $k (keys %ParamDescription) {
      if ( $s =~ /$k/i ) {
        $key = $k;
        last;
      }
    }
    AssignParam($key,$s) if defined($key);
  }
  close PARFILE;

  # Catch errors in the parameter file.
  if ( ($MinAirWn <= 0) and ($MaxAirWn <= 0) ) {
    $MinAirWn  =  1e15;
  }
  if ($IntLastCol < $IntFirstCol) {
    Error("Error in value of last column of line intensity in par. file. It must be >= the 1st column of line intensity.");
  }
  $widthInt = $IntLastCol - $IntFirstCol + 1;
  if ($WLLastCol < $WLFirstCol) {
    Error("Error in value of last column of line wavelength in par. file. It must be >= the 1st column of line wavelength.");
  }
  if ($WuncLastCol < $WuncFirstCol) {
    Error("Error in value of last column of line uncertainty in par. file. It must be >= the 1st column of line uncertainty.");
  }
  if ($WuncLastCol < $WuncFirstCol) {
    Error("Error in value of last column of line weight in par. file. It must be >= the 1st column of line weight.");
  }
  if ($LLLastCol < $LLFirstCol) {
    Error("Error in value of last column of lower level label in par. file. It must be >= the 1st column of lower level label.");
  }
  $LNameWidth = $LLLastCol - $LLFirstCol + 1;

  # Assign default values to missing parameters
  my $missing = '';
  foreach my $k (keys %ParamDescription) {
    my ($var_name,$type,$valid_range,$default_value,$is_defined) = @{$ParamDescription{$k}};
    if ( !$is_defined ) {
      if ( !defined($default_value) ) {
        $missing .= "$k mandatory parameter missing in the par. file.\n";
      } else {
        eval("$var_name=$default_value");
      }
    }
  }
  Error($missing) if $missing;
}

sub PrintMatrix($$) {
  # Write the diagnostics file LOPT.OUT in case of degenerate matrix.
  my ($G,$NA) = @_;
  my $NCol = 5;
  my ($i,$j,$i1,$j1,$k,$jmin,$jmax);

  print 'Printing matrix to the LOPT.OUT file ...';
  open(FOut,'>LOPT.OUT') or
    Error('Error opening output file.');
  print FOut "Matrix A(i,j) (dimension=$NA) :\n", MatrixToString($A,"\t") or
        Error('Error writing to output file.');
  print FOut "Vector E(j) (dimension=$NA :\n";
  for (my $k = 1; $k <= ($NA-1) / $NCol + 1; $k++) {
    my $jmin = ($k-1)*$NCol + 1;
    my $jmax = $jmin+$NCol-1;
    $jmax = $NA if ($jmax>$NA);
    for (my $j = $jmin; $j <= $jmax; $j++) {
      print FOut sprintf("%10d    ",$j),$E[$j-1]+0,"\n";
    }
  }
  close(FOut) or Error('Error closing output file.');
  print "Done.\n";
}

############################################################################
sub Substr($$;$) {  #09/29/2009 9:37AM
############################################################################
  # Return a substring of a string (first parameter),
  # starting from the given position (second parameter),
  # with a given length (optional third parameter).
  # If length of the substring is omitted, return everything to the end of
  # the string.
  my ($s,$start,$len) = @_;
  my $strlen = length($s);
  $len = $strlen - $start unless $len;
  return '' if ( ($strlen < 0) || ($len < $0) || ($start > $strlen) );
  $len = $strlen - $start if $start + $len > $strlen;
  return substr($s,$start,$len);
} ##Substr($$;$)

############################################################################
sub leastSignificantFigureValue($$$) {   #11/12/2009 7:33AM
############################################################################
  # Return the value of the least significant digit of a properly rounded
  # floating-point number x having uncertainty unc,
  # where rounding is made using the rule-of-m rounding with
  # m = given round-off threshold (third parameter).
  my ($x,$unc,$round_thresh) = @_;
  $x = sprintf("%.20f", $x);
  $x =~ s/^\s+|\s+$//g;
  my $s_unc = sprintf("%.10e",$unc); # This rounds $unc to 11 significant digits
  $s_unc =~ /^([^e]+)e(.+)$/i;
  my $unc1 = $1*100;     # between 100.0 and 999.9
  my $pl =  $unc1 > $round_thresh*100 ? 3 : $unc1 > $round_thresh*10 ? 2 : $unc1 > $round_thresh ? 1 : 0;
  my $power = $2 + $pl - 2;
  my $s_d = "1e$power";
  return $s_d;
} ##leastSignificantFigureValue($$$)

my $_log10 = 1/log(10);
############################################################################
sub get_power($$) {  #11/13/2009 9:04AM
############################################################################
  # Return the exponent in the decimal scientific format.
  # This function is used in the round() function below.
  my ($unc, $round_thresh) = @_;
  my $s_unc = sprintf("%.10e",$unc); # This rounds $unc to 11 significant digits
  $s_unc =~ /^([^e]+)e(.+)$/i;
  my $pl = int(log($1/$round_thresh)*$_log10);
  my $power = $2 + $pl;
  return $power;
} ##get_power($$)

############################################################################
sub round($$$$) {   #09/30/2009 3:37PM
# The first and second arguments are the value to be rounded and its uncertainty (also to be rounded).
# The third argument is the precision (a.k.a. round-off threshold), i.e. the maximum uncertainty in units of the least significant digit.
# The fourth argument is a flag denoting whether or not we want to add the actual rounding error to the uncertainty value.
# Returns the rounded value, the rounded uncertainty, and the value of the least significant digit (all as character strings).
############################################################################
  my ($x,$unc,$round_thresh,$add_rounding_error) = @_;
  # $round_thresh must be greater than 2 - in the units of the least significant digit

  my ($s_x, $diff);

  $unc = abs($unc);

  if ( $unc < $too_small ) {
    return ($x,'_');
  }

  $x = sprintf("%.20f", $x);
  $x =~ s/^\s+|\s+$//g;

  my $power = get_power($unc,$round_thresh);
  my ($lsd,$s_unc);
  if ( $power >= 0 ) {
    my ($factor1, $factor2) = ("1e$power","1e-$power");
    $s_x   = sprintf("%.0f",$x*$factor2)*$factor1;
    if ( $add_rounding_error ) {
      $diff = $x - $s_x;
      $unc  = sqrt($unc*$unc + $diff*$diff);
      $power = get_power($unc,$round_thresh);
      ($factor1, $factor2) = ("1e$power","1e-$power");
      $s_x   = sprintf("%.0f",$x*$factor2)*$factor1;
    }
    $s_unc = sprintf("%.0f",$unc*$factor2)*$factor1;
    $lsd = "1e$power";
  } else {
    $lsd = "1e$power";
    $power = -$power;

    my $dot_pos = index($x,'.');
    my $s_x_int = ($dot_pos >= 0) ? substr($x,0,$dot_pos) : $x;
    my $s_x_decimals = ($dot_pos >= 0) ? substr($x,$dot_pos) : '';
    $s_x_decimals = sprintf("%.${power}f",$s_x_decimals);
    # Bug fix July 2 2017 by M. Ya. Vokhmentsev
    #$s_x_int++ if substr($s_x_decimals,0,1) > 0;
    $s_x_int += ($x >= 0) ? 1 : -1 if substr($s_x_decimals, 0, 1) > 0;
    # End fix July 2 2017 b
    $s_x   = $s_x_int . substr($s_x_decimals,1);

    if ( $add_rounding_error ) {
      $diff = $x - $s_x;
      # Combine the initial uncertainty with the rounding error.
      # The estimate below is close to reality when $round_thresh is
      # greater than or equal to 5. Otherwise, since the rounding errors
      # have uniform distribution instead of normal (Gaussian) one,
      # the rounding errors result in significant distortion
      # of distribution of the results and, in addition,
      # produce systematic errors depending on the last digits of the
      # quantity-to-be-rounded.
      $unc  = sqrt($unc*$unc + $diff*$diff);
      $power = get_power($unc,$round_thresh);
      $lsd = "1e$power";
      if ( $power >= 0 ) {
        my ($factor1, $factor2) = ("1e$power","1e-$power");
        $s_x   = sprintf("%.0f",$x*$factor2)*$factor1;
        $s_unc = sprintf("%.0f",$unc*$factor2)*$factor1;
      } else  {
        $power = -$power;
        $dot_pos = index($x,'.');
        $s_x_int = ($dot_pos >= 0) ? substr($x,0,$dot_pos) : $x;
        $s_x_decimals = ($dot_pos >= 0) ? substr($x,$dot_pos) : '';
        $s_x_decimals = sprintf("%.${power}f",$s_x_decimals);
        # Bug fix July 2 2017 by M. Ya. Vokhmentsev
        #$s_x_int++ if substr($s_x_decimals,0,1) > 0;
        $s_x_int += ($x >= 0) ? 1 : -1 if substr($s_x_decimals, 0, 1) > 0;
        # End fix July 2 2017 b
        $s_x   = $s_x_int . substr($s_x_decimals,1);
        $s_unc = sprintf("%.${power}f",$unc);
      }
    } else {
      $s_unc = sprintf("%.${power}f",$unc);
    }
  }
  $s_x =~ s/^[ -]+// if $s_x == 0;
  return ($s_x,$s_unc,$lsd);
} ##round($$$$)

############################################################################
sub num_dec_places($) {  #10/13/2009 12:12PM
############################################################################
  # Return the number of places after the decimal point in a floating-point number
  my $X = shift;
  my $dot_index = index($X,'.');
  my $decimalPlaces = ($dot_index >=0) ? length($X) - $dot_index - 1 : 0;
  return $decimalPlaces;
} ##num_dec_places($)

sub PrintRound($$$) {
  # Return a formatted decimal-point number, aligned by the decimal point
  # within a given width, left- and right-padded by spaces.
  my ($Xr,$width,$WidthAfterPoint) = @_;

  return $Xr if ($TabDelimitedOutput);
  return NotAvail($width,$WidthAfterPoint) if ($Xr eq '_');

  my $places = num_dec_places($Xr);

  $Xr = LeftPad($Xr,' ',$width);
  my $j1 = $WidthAfterPoint - $places;
  $j1++ unless ($places || ($WidthAfterPoint <= 0));
  $j1 = 0 if $j1 < 0;
  for (my $j2 = 1; $j2 <= $j1; $j2++) {
    # Move space from the beginning of the string to its end
    $Xr = substr($Xr,1) . ' ' if (substr($Xr,0,1) eq ' ');
  }
  return $Xr;
}

sub FindNumAssign {
  # Calculate the number of assigned transitions for each observed spectral
  # line and the weights of transitions in multiply-assigned lines.
  # Divide the measurement uncertainty of the transition by the square
  # root of its weight in the blend. This is equivalent to
  # multiplying the transition weight in the least-squares optimization by
  # the weight in the blend.
  foreach my $wl ( keys %LinesByW ) {
    my ($TotWeight, $MaxWeight, $NumAssigned) = (0,0,0);
    if ( $wl eq '328.090' ) {
      $wl = $wl;
    }
    foreach my $PT (@{$LinesByW{$wl}}) {
      if ($PT->{'Weight'} > $too_small) {
        $NumAssigned++;
        $TotWeight += $PT->{'Weight'};
        $MaxWeight = $PT->{'Weight'} if ($MaxWeight < $PT->{'Weight'});
      }
    }
    foreach my $PT (@{$LinesByW{$wl}}) {
      if ($PT->{'Weight'} > $too_small) {
        $MaxWeight /= $TotWeight;
        $PT->{'Weight'} = sqrt($PT->{'Weight'}/$TotWeight);
        # Divide the uncertainty by square root of weight in the blend.
        $PT->{'WnUnc'} /= $PT->{'Weight'};
        $PT->{'NumAssign'} = $NumAssigned;
      }
    }
  }
  # Destroy the sorted lines hash to free memory.
  %LinesByW = undef;
}

sub ReadLines {
      # Read transitions data from the transitions-input file

  open(Tfil,"<$TransFilename") or
    Error("Unable to open transitions input file $TransFilename");

  $LineNo = 0;
  $LastWunit = 'A'; # Default line unit is Angstrom
  $LastWLUnc = '';
  my $WLfactor = 1;
  $MinUnc = 1e100;
  $MaxCGNo = 0;
  my $Res = 0;
  my $S;
  while (($S = <Tfil>) && defined($S)) {
    chomp($S);
    my $PT = InitTrans();
    my $PL1 = InitLevel();
    my $PL2 = InitLevel();
    $PT->{'FromLevel'} = $PL2;
    $PT->{'ToLevel'} = $PL1;
    $LineNo++;
     #  Read line flags
    my $S1 = Substr($S,$FlgFirstCol-1,5);
    $S1 =~ s/^\s+|\s+$//g; # Trim spaces
    if ( $S1 =~ /^([^\d]*)(\d+)([^\d]*)$/ ) {
      my ($digits, $non_digits) = ($2,$1.$3);
      $PT->{'CGroupNo'} = $digits;
      $MaxCGNo = $digits if ( $digits > $MaxCGNo );
      $S1 = $non_digits;
    }
    $Res = 0;

    my $flg = 0;
    $flg += $Fair if ($S1 =~ /A/i);
    $flg += $FBlend if ($S1 =~ /B/i);
    $flg += $FQuest if ($S1 =~ /Q/i);
    $flg += $FMask if ($S1 =~ /M/i);
    $flg += $FPredict if ($S1 =~ /P/i);
    $PT->{'Flags'} = $flg;

      # Read wavelength units
    my $WunitStr = Trim(Substr($S,$WunitFirstCol-1,4));
    foreach my $key (keys %WLfactors) {
      $WunitStr =~ s/^$key$/$key/i;
    }

    if (!$WunitStr) {
      $WunitStr = $LastWunit; # Angstroms by default
    } else {
      if (!defined $WLfactors{$WunitStr}) {
        Error("Error in wavelength units in transition-input file, line No. $LineNo");
      }
    }
    $PT->{'Wunit'} = $WunitStr;
    $LastWunit = $WunitStr;
    $WLfactor = $WLfactors{$WunitStr};

     # Read the input quantity (wavelength or wavenumber or frequency or energy interval) value
     # and convert to cm-1
    $S1 = Trim(Substr($S,$WLFirstCol-1,$WLLastCol-$WLFirstCol+1));
    $PT->{'Wobs'} = $S1;
    my $point_pos = index($S1,'.');
    my $decimals = ($point_pos >= 0 && $point_pos < length($S1)) ? substr($S1,$point_pos+1) : '';
    my $num_decimals = length($decimals);
    my $lambda = $S1 + 0;
    $lambda = 0.0 unless ValidNumber($S1);
    $Res = -1 if ((abs($lambda) < $too_small) and ($WLfactor>0) or !ValidNumber($S1));
    my $precision = ($Res) ? 0 : "1e-$num_decimals" / ($lambda ==0 ? 1e-27 : $lambda);

    if ( !($PT->{'Flags'} & $FMP) || ($lambda > $too_small) ) {
      if ($Res) {
        Error("Error in wavelength in transition-input file, line No. $LineNo");
      }

          # Convert lambda to angstroms
      if ($WLfactor >= 0) {
        $lambda *= $WLfactor;
        if ( $PT->{'Flags'} & $Fair ) {
          $PT->{'Eobs'} = 1e8/Lvac(abs($lambda));
          $PT->{'Eobs'} = -$PT->{'Eobs'} if ( $lambda < 0 );
        } else {
          $PT->{'Eobs'} = 1e8/$lambda;
          if ( AirWavenumber($PT->{'Eobs'}) ) {
            $PT->{'Flags'} |= $Fair;
            $PT->{'Eobs'} = 1e8/Lvac(abs($lambda));
            $PT->{'Eobs'} = -$PT->{'Eobs'} if ($lambda < 0);
          }
        }
      } else {  # Wavenumbers and frequencies are always in vacuo
          # Make sure the $Fair flag is cleared
        $PT->{'Flags'} |= $Fair;
        $PT->{'Flags'} -= $Fair;

        $PT->{'Eobs'} = -1e8*$lambda/$WLfactor;
        if (abs($lambda) < $too_small) {
          $lambda = 0.0;
        } else {
          $lambda = -$WLfactor/$lambda;
        }
      }

    }

     # Read Intensity (empty entries can be replaced by "/" or "_"
    $S1 = Substr($S,$IntFirstCol-1,$IntLastCol-$IntFirstCol+1);
    $S1 =~ s/\s+$//g; # Trim only the trailing spaces.
    $PT->{'Iobs'} = $S1;

     # Read uncertainty of the input quantity and convert to cm-1
    $S1 = Trim(Substr($S,$WuncFirstCol-1,$WuncLastCol-$WuncFirstCol+1));
    unless ( $S1 ) {
      $S1 = $LastWLUnc;
    }
    $LastWLUnc = $S1;
    $PT->{'WnUnc'} = $S1;
    if ( $PT->{'Flags'} & $FMP ) {
      $PT->{'WnUnc'} = $S1;
    }
    if (!ValidNumber($S1) or ($PT->{'WnUnc'} < $too_small)) {
      if (!($PT->{'Flags'} & $FMP)) {
        Error("Error in wavelength uncertainty in transition-input file. Line No. $LineNo");
      } else {
        $PT->{'WnUnc'} = $PT->{'Eobs'} * $precision;
      }
    }
    if (!($PT->{'Flags'} & $FMP) or ($lambda > $too_small)) {
      if ($lambda>1e-10) {
        if ($WLfactor>0) {
          $PT->{'WnUnc'} = 1.e8*$PT->{'WnUnc'}*$WLfactor/sqr($lambda);
        } else {
          $PT->{'WnUnc'} = -1.e8*$PT->{'WnUnc'}/$WLfactor;
        }
      }
    }
    $MinUnc = $PT->{'WnUnc'} if ( ($PT->{'WnUnc'} < $MinUnc) && !($PT->{'Flags'} & $FMP) );

     # Read line weight - for multiply assighned lines
    $S1 = Trim(Substr($S,$WeightFirstCol-1,$WeightLastCol-$WeightFirstCol+1));
    if ($PT->{'Flags'} & $FMP) {
      $PT->{'Weight'} = 0;
    } elsif ($S1) {
      $PT->{'Weight'} = $S1 + 0;
      if ((!ValidNumber($S1)) or ($PT->{'Weight'} < $too_small)) {
        Error("Error in line weight in transition-input file. Line No. $LineNo");
      }
    }

     # Read lower level label
    $PL1->{'Name'} = Trim(Substr($S,$LLFirstCol-1,$LNameWidth));

     # Read upper level label
    $PL2->{'Name'} = Trim(Substr($S,$ULFirstCol-1,$LNameWidth));

     # Insert the levels in the level group lists. Merge groups
     # connected by this transition.
    InsertLevs($PT);
     # Add the transition to the Lines array
    push(@Lines,$PT);
     # Add the line to the hash sorted by wavelength
    my $key = $PT->{'Wobs'} . $PT->{'Wunit'};
    $LinesByW{$key} = [] unless defined $LinesByW{$key};
    push(@{$LinesByW{$key}}, $PT);
     # Add the line to the hash sorted by lower and upper levels
    $key = $PT->{'ToLevel'}->{'Name'} . ',,' . $PT->{'FromLevel'}->{'Name'};
    $LinesByLevels{$key} = $LineNo;
     # Add the transition to the corresponding arrays of the lower and upper levels
    push(@{$PT->{'ToLevel'}->{'To'}},$PT) unless ($PT->{'Flags'} & ($FMP));
    push(@{$PT->{'FromLevel'}->{'From'}},$PT) unless ($PT->{'Flags'} & ($FMP));
  }
  close(Tfil);
}

sub FillCGroups {
  # Initialize the arrays associated with correlated line groups, if any.
  return unless $MaxCGNo;

  @CGroups = ();
  for (my $i = 1; $i <= $MaxCGNo; $i++) {
    push(@CGroups, []);
  }
  $CUnc = [];
  for (my $i = 1; $i <= $MaxCGNo; $i++) {
    push(@{$CUnc}, 0);
  }

  foreach my $PT (@Lines) {
       # For each line, initialize the output array of uncertainties
    my $C = [];
    for (my $i = 1; $i <= $MaxCGNo; $i++) {
      push(@{$C}, 0);
    }
    $PT->{'C'} = $C;

    my $j = $PT->{'CGroupNo'};
    if ($j && !($PT->{'Flags'} & $FMP)) {
      my $PCG = $CGroups[$j-1];
      push(@{$PCG}, $PT);
      # set $CUnc->[$j-1] to minimum measured line uncertainty in the line group
      if ( !$CUnc->[$j-1] or (abs($PT->{'WnUnc'}) < $CUnc->[$j-1]) ) {
        $CUnc->[$j-1] = abs($PT->{'WnUnc'});
      }
    }
  }
    # For each level, initialize the output array of uncertainties
  foreach my $PL (@Levels) {
    my $C = [];
    for (my $i = 1; $i <= $MaxCGNo; $i++) {
      push(@{$C}, 0);
    }
    $PL->{'C'} = $C;
  }
}

sub InsertVirtLine($$$) {
  my ($Level,$RefLevel,$D)= @_;
  # Insert virtual lines between fixed levels and the ground state
  # or between pairs of fixed levels if their uncertainty is given
  # relative to each other
  #

  return if ($Level eq $RefLevel);

  my $PT = InitTrans();
  $PT->{'Flags'} = $FVirtual + (($D eq '') && defined($RefLevel) ? $FPredict : 0);
  $PT->{'FromLevel'} = $Level;
  $PT->{'ToLevel'} = $RefLevel;
  if (!defined $RefLevel) {
    $Level->{'D'} = $D;       # uncert. relative to ground
  } else {
    $PT->{'Eobs'} = abs($Level->{'E'}-$RefLevel->{'E'});
    if ($Level->{'E'} < $RefLevel->{'E'}) { # invert levels to make Eobs positive
      $PT->{'FromLevel'} = $RefLevel;
      $PT->{'ToLevel'} = $Level;
    }
  }
  $PT->{'WnUnc'} = $D;

  if (defined($RefLevel) and (!$DivideGroups || ($RefLevel->{'D'} != 0))) {
    InsertLevs($PT);
  }
  push(@Lines,$PT);
  push(@VirtLines,$PT) unless ($PT->{'Flags'} & $FPredict);
}

sub EvalFixedLevUncert($$$$$) {
  # Parse the uncertainty specification of a fixed level,
  # as read from the input fixed-levels file.
  my($PL,$S,$D,$PLref,$LNo) = @_;
  my ($j,$Res,$S2,$FromGround,$LabelDef,$Ref);
  my ($Ev,$must_insert_virt_line) = (0,0);

  $S =~ s/^\s+|\s+$//g; # Trim
  return ($S,$D,undef,0) if ($S eq '');

  $LabelDef =  ($S =~ /^\[/);
  $S =~ s/^[^0-9. ]+//g; # Delete any leading non-digits
  return ($S,$D,undef,0) if ($S eq '');

  my $S2 = '';
  if ( $S =~ /^([0-9.]+)([^0-9. ]|$)/ ) {
    $S2 = $1;
    $S =~ s/^$S2//;
  }
  $S =~ s/^\s+|\s+$//g; # Trim
  $FromGround = ( !$S or ($S !~ /^\(/) );
  $S2 =~ s/^\s+|\s+$//g; # Trim
  if ( $S and ($S !~ /^\[/) ) {
    $S =~ s/^.//;  # Delete the leading symbol, which can only be '(' or ']'
  }
  $S =~ s/^\s+|\s+$//g; # Trim
  if ($LabelDef) {
    $Ref = $S2 + 0;
    $PLref = undef;
    foreach my $PL1 (@FixLevs) {
      if ( $PL1->{'Fixed'} == $Ref ) {
        ($PLref = $PL1) && last;
      }
    }
    if (defined($PLref)) {
      Error("Duplicate fixed level label, line No. $LNo");
    }
    $PL->{'Fixed'} = $Ref;
    return ($S,$D,$PLref,0);
  }
  $D = $S2 + 0;
  if ( ($D < 0) || ($S2 !~ /^\d*\.{0,1}\d*/) ) {
    Error("Error in fixed level uncertainty, line No. $LNo");
  }
  $MinUnc = $D if ( ($D < $MinUnc) && ($D > 0) );
  $PLref = undef;
  if (!$FromGround) { # Get the reference level
    $S2 = '';
    if ( $S =~ /^(\d+)([^\d]|$)/ ) {
      $S2 = $1;
      $S =~ s/^$S2//;
      $S =~ s/^.//; # Delete the leading symbol, which should be ')'
    }
    $S2 =~ s/^\s+|\s+$//g; # Trim
    $S =~ s/^\s+|\s+$//g; # Trim
    $Ref = $S2 + 0;
    $PLref = undef;
    foreach my $PL (@FixLevs) {
      if ( $PL->{'Fixed'} == $Ref ) {
        ($PLref = $PL) && last;
      }
    }
    if (!defined($PLref)) {
      Error("Undefined reference level label for fixed level, line No. $LNo");
    }
    $must_insert_virt_line = 1 unless $PL->{'Free'};
  } else {
    $PL->{'D'} = $D;
    $must_insert_virt_line = 1 unless $PL->{'Free'};
  }
  return ($S,$D,$PLref,$must_insert_virt_line);
}

sub ReplaceGround($) {
  # For a virtual line created in the ReadFixLevs procedure,
  # if the lower level ('ToLevel') is not defined, set the lower level to
  # $Ground. Insert the upper and lower levels in the level groups and
  # combine the groups if necessary.
  my $PT = shift;
  if (($PT->{'Flags'} & $FVirtual) and !$PT->{'ToLevel'}) {
    $PT->{'ToLevel'} = $Ground;
    if (!$DivideGroups) {
      InsertLevs($PT);
    }
    $PT->{'Eobs'} = $PT->{'FromLevel'}->{'E'} - $Ground->{'E'};
  }
}

sub ReadFixLevs {
  # Read the input fixed-levels file
  my ($S,$S1,$FLName,$extra,$LNo,$NL,$E,$U,$G,$PL,$PLref,$PT);

  open(FixLF,"<$FixLevFilename") or
    Error("Unable to open fixed levels file $FixLevFilename");

  $Ground = undef;
  @FixLevs = ();

  $LNo = 0;
  while (($S = <FixLF>) && defined($S)) {
    $LNo++;
    chomp($S);
    $extra = '';
    $S =~ s/;.*$//g;  # delete comments (to the right of a semicolon)
    $S =~ s/\s+$//g; # delete trailing spaces
    next unless $S;   # Skip blank lines and lines starting with semicolons

    if ( $S =~ /^(.{$LNameWidth})\s*([0-9.+-e]+)\s+([0-9.+-e\[\]\(\)]+)(\s+[^;]+)*/i ) {
      ($FLName,$E,$S1,$extra) = ($1,$2,$3,$4);
    } else {
      Error("Error in fixed level file format, line No. $LNo");
    }

    $FLName = Trim($FLName);

    # Verify the fixed level energy value
    if ( $E =~ /^[+-]{0,1}[0-9]+[.]{0,1}[0-9]*(e[+-]{0,1}\d+){0,1}$/ ) {
      $E = $E + 0;
    } else {
      Error("Error in fixed level energy, line No. $LNo");
    }

    $extra =~ s/^\s+|\s+$//g; # delete leading and trailing spaces
    $extra = lc($extra);

    # Locate the fixed level in the groups read from the transitions file
    $PL = undef;
    for (my $i = 1; $i <= $#{$Groups} + 1; $i++) {
      my $G = $Groups->[$i-1];
      foreach my $PLx (@{$G}) {
        if ($PLx->{'Name'} eq $FLName) {
          $PL = $PLx;
          $PL->{'E'} = $E;
          $PL->{'GNo'} = $i;
          last;
        }
      }
      last if defined($PL);
    }

    if (!defined($PL)) {
      print "\nWarning: " .
       ( $extra =~ /base/ ? 'Base level' : 'Level to fix') .
        " [$FLName] is not found in the transitions list;" .
        " fixed-levels file line No. $LNo\n";

      my ($G1,$k1) = defined($Levels{$PL->{'Name'}}) ? @{$Levels{$PL->{'Name'}}} : (undef,-1);
      if ( defined($G1) ) {
        $PL = $G1->[$k1];
      } else {
        $PL = InitLevel();
        $PL->{'Name'} = $FLName;
        $PL->{'E'} = $E;
        my $G2 = [$PL];
        push(@{$Groups},$G2);
        $PL->{'GNo'} = $#{$Groups} + 1;
        $Levels{$PL->{'Name'}} = [$G2,$#{$G2}];
        push(@Levels,$PL);
      }
    }

    #if ( $PL->{'Name'} eq '74') {
    #  $PL = $PL;
    #}

    if ( ($extra =~ /ground/i) || (($E == 0) && ($extra !~ /base|free/i)) ) {
      if ( $Ground ) {
        Error("Error in fixed-levels line No. $LNo: duplicate ground level defined")
      }
      $Ground = $PL;
    }
    if ( ($extra =~ /base/i) && ($PL ne $Ground) ) {
      push(@UncBases,$PL);
    }
    if ( $extra =~ /free/i ) {
      $PL->{'Free'} = 1;
    }

    if ( $E || ($PL eq $Ground) ) {
      push(@FixLevs,$PL);

      # Interpret fixed level uncertainty:
      # Insert virtual lines according to specified uncertainty options
      my $must_insert_virt_line = 0;
      do {
        ($S1,$U,$PLref,$must_insert_virt_line) = &EvalFixedLevUncert($PL,$S1,$U,$PLref,$LNo);
        if ($must_insert_virt_line) {
          InsertVirtLine($PL,$PLref,$U);
        }
      } while ($S1 ne '');
    }
  }
  close(FixLF);

  if ( !defined($Ground) ) {
    # Locate and fix the ground level
    foreach my $FL (@FixLevs) {
      if (!defined($Ground)
        || ($FL->{'D'} < $Ground->{'D'})
        || (($FL->{'D'} == $Ground->{'D'}) && ($FL->{'E'} < $Ground->{'E'}))
      ) {
        $Ground = $FL;
      }
    }
  }

  if (defined($Ground)) {
    $Ground->{'D'} = 0;
  } else {
    Error("Could not find the ground level among the fixed levels.");
  }

  # Set Fixed flags for fixed levels
  foreach my $FL (@FixLevs) {
    if ( ($FL->{'D'} > 99999) || $FL->{'Free'} ) {
      $FL->{'Fixed'} = 0;
    } else {
      $FL->{'Fixed'} = 1;
    }
    $FL->{'D'} = 0;
  }

  for (my $ib = 0; $ib <= $#UncBases; $ib++) {
    my $BL = $UncBases[$ib];
    my $GN = $BL->{'GNo'};
    my $G = $Groups->[$GN-1];
    foreach my $Lev (@{$G}) {
      $Lev->{'UncRelBases'} = [] unless defined($Lev->{'UncRelBases'});
      $Lev->{'UncRelBases'}->[$ib] = 0;
      InsertVirtLine($BL,$Lev,'');
    }
  }

  # Make sure at least one level is fixed in each group
  my $GNo = 0;
  my $UnboundDetected = 0;
  foreach my $G (@{$Groups}) {
    $GNo++;
    next unless (defined($G) && $G->[0]);
    my $PL = undef;
    foreach my $GL (@{$G}) {
      if ($GL->{'Fixed'}) {
        $PL = $GL;
        last;
      }
    }
    if (!defined($PL)) {  # Group is unbound; fix the first level
      $UnboundDetected++;
      $DivideGroups = 1;
      $PL = $G->[0];
      $PL->{'Fixed'} = 2;
      $PL->{'E'} = 0;
      $PL->{'GNo'} = $GNo;
      if (!defined($Ground)) {
        $Ground = $PL;
        $PL->{'D'} = 0;
      } else {
        $PL->{'D'} = 0.0001;
      }
      push(@FixLevs, $PL);
      InsertVirtLine($PL,$Ground,$PL->{'D'});
    }
  }

  if ( $UnboundDetected) {
    print "\nWarning: $UnboundDetected unbound groups detected. " .
    ($DivideGroups ? '' : 'Switching to option "Divide levels into independent groups = Y".' ."\n") .
    "One level in each unbound group is fixed at zero energy.\n";
  }

  foreach my $PT (@Lines) {
    &ReplaceGround($PT);
  }

  # Remove transitions having the same upper and lower levels
  my $repeat = 1;
  while ($repeat) {
    $repeat = 0;
    for (my $i = 0; $i <= $#Lines; $i++) {
      my $PT = $Lines[$i];
      if ($PT->{'FromLevel'} eq $PT->{'ToLevel'}) {
        splice(@Lines,$i,1);
        $repeat = 1;
        last;
      }
    }
  }

    # Reorder each group so that all fixed levels are in the end
    # of the group
  for (my $i=1; $i <= $#{$Groups} + 1; $i++) {
    my $G = $Groups->[$i-1];
    foreach my $PL (@{$G}) {
      $PL->{'GNo'} = $i;
    }
    my $j = 0;
    my $NL = $#{$G} + 1;
    while ($j < $NL) {
      $j++;
      my $PL = $G->[$j-1];
      if ($PL->{'Fixed'} && ($j < $NL)) {
        my $PL1 = $G->[$NL-1];
        while (($NL > $j) && $PL1->{'Fixed'}) {
          $NL--;
          $PL1 = $G->[$NL-1];
        }
        if (($NL>$j) && !$PL1->{'Fixed'}) {
          # Switch $PL and $PL1
          $G->[$j-1]  = $PL1;
          $G->[$NL-1] = $PL;
        }
      }
    }
  }

  #$FixedLevUncert = $MinUnc*$FixWeightFactor;
  #$FixedLevWeight = 1/$FixedLevUncert**2; # Weight used to effectively fix the "fixed" levels.
}

############################################################################
sub IndexOf($$) {  #12/03/2006 11:35AM A.Kramida
############################################################################
  # Return the index of a given object (first parameter)
  # in a given array (second parameter).
  # Return -1 if not found.
  my ($P, $A) = @_;
  my $ind = -1;
  foreach (my $i = 0; $i <= $#{$A}; $i++) {
    if ( $P eq $A->[$i] ) {
      $ind = $i;
      last;
    }
  }
  return $ind;
} ##IndexOf($$)

sub InsertFixedLev($$) {
  my ($PL,$G) = @_;
  return unless ($PL->{'Fixed'} == 1);
  # Find this level in group $G
  foreach my $L (@{$G}) {
    if ( $PL eq $L) {
      return;
    }
  }
  # Level was not found, so we insert it in group $G
  push(@{$G}, $PL);
}

############################################################################
sub GetLineWeight($$) {   #11/24/2009 2:42PM
############################################################################
  # Return the weight of the given line (first parameter) in the
  # level-optimization problem. The weight of virtual lines is different
  # for the two level optimizations done in Pass==1 and Pass==2 (second parameter).
  my ($PT,$Pass) = @_;
  if (!($PT->{'Flags'} & $FMP) && ($PT->{'WnUnc'} < $too_small)) {
    Error('Error: zero line uncertainty.');
  }
  my $FF = ($PT->{'Flags'} & ($FMP | $FExcluded)) ? 0 : 1/sqr($PT->{'WnUnc'});

  if (($PT->{'Flags'} & $FVirtual) && ($Pass==1) && ($PT->{'WnUnc'} <= 99999) ) {
    #$FF = $FixedLevWeight;    # This should really fix the level
    $FF /= $FixWeightFactor;
  }
  return $FF;
} ##getLineWeight($$)

my $num_subst = 0;

############################################################################
sub FindSubstitutions($$$) {  #11/24/2009 2:10PM
############################################################################
  # Parameters: Levels group $G, stage of level optimization $Pass,
  # number of substitutions from the previous stage $num_subst.
  # Returns number of substitutions.
  #
  # Variable substitutions may be necessary in order to reduce effects of
  # numerical cancellation which may lead to large numerical errors in the
  # construction of the matrix elements as well as in the matrix inversion
  # procedure. The substitution consists in replacing the unknown energy
  # level variables by smaller and more precisely defined differnces between
  # them. Read more about the algorithm in the CPC article.
  # The substitutions are defined in the lists attached to each substituted
  # level.
  my ($G,$Pass,$num_subst) = @_;
  foreach my $L (@{$G}) {
    #next if ($L eq $Ground);
    $L->{'Subst'} = undef;
  }
  return 0 unless $use_substitution;
  return $num_subst if $Pass > 1;

  # A. Kramida Aug 24 2017: Lines related to @virt and @virtThisGroup
  # are a bug fix intended to make variable substitution work with
  # various flavors of fixed levels and with division into independent groups
  my @virtThisGroup = ();
  foreach my $PT (@VirtLines) {
    next if ( $PT->{'FromLevel'} eq $PT->{'ToLevel'} );
    next if ( $PT->{'FromLevel'}->{'CurGNo'} != $Ground->{'CurGNo'} );
    next if ( $PT->{'ToLevel'}->{'CurGNo'} != $Ground->{'CurGNo'} );
    push (@virtThisGroup,$PT);
  }
  foreach my $L (@{$G}) {
    #next if ($L->{'Fixed'}); # Commented as part of the bug fix of Aug 24 2017
    next if ($L eq $Ground);
    my @virt  = ();
    foreach my $PT (@virtThisGroup) {
      push (@virt,$PT) if (($PT->{'FromLevel'} eq $L) || ($PT->{'ToLevel'} eq $L));
    }
    # Find max and min weights of transitions to or from each level
    my $min_weight = 1e100;
    my $max_weight = 0;
    my $T_max = undef;
    #foreach my $T (@{$L->{'To'}},@{$L->{'From'}}) {  # Commented as part of the bug fix of Aug 24 2017
    foreach my $T (@{$L->{'To'}},@{$L->{'From'}},@virt) {
      next if ($T->{'Flags'} & $FExcluded);
      my $weight = GetLineWeight($T,$Pass);
      $min_weight = $weight if ($weight < $min_weight);
      if ($weight > $max_weight) {
        $max_weight = $weight;
        $L->{'max_weight'} = [$T,$max_weight];
        $T_max = $T;
      }
    }
    #next if $T_max->{'FromLevel'}->{'Fixed'}; # Commented as part of the bug fix of Aug 24 2017
    #next if $T_max->{'ToLevel'}->{'Fixed'};   # Commented as part of the bug fix of Aug 24 2017
    next if $T_max->{'ToLevel'} eq $Ground;
    next if $T_max->{'FromLevel'} eq $Ground;

    # If max and min weights differ too much, the substitutions are necessary.
    my $threshold = 1e5;
    if ( $max_weight/$min_weight > $threshold ) {
      my $L_from = $T_max->{'FromLevel'};
      my $L_to = $T_max->{'ToLevel'};
      if ( (!defined($L_from->{'Subst'}) || ($L_from->{'Subst'}->[1] < $max_weight)) &&
        (!defined($L_from->{'max_weight'}) ||
          #(!$L_from->{'max_weight'}->[0]->{'ToLevel'}->{'Fixed'}) ||  # Commented as part of the bug fix of Aug 24 2017
          ($L_from->{'max_weight'}->[1] < $max_weight)
        )
      ) {
        $num_subst++ if ( !defined($L_from->{'Subst'}) );
        $L_from->{'Subst'} = [$L_to,$max_weight];
      }
      foreach my $T (@{$L->{'To'}}) {
        next if ($T->{'Flags'} & $FExcluded);
        my $weight = GetLineWeight($T,$Pass);
        if ( $weight/$min_weight > $threshold) {
          $num_subst++ if ( !defined($T->{'FromLevel'}->{'Subst'}) );
          if ( (!defined($T->{'FromLevel'}->{'Subst'}) || ($T->{'FromLevel'}->{'Subst'}->[1] < $weight)) &&
            (!defined($T->{'FromLevel'}->{'max_weight'}) ||
              #(!$T->{'FromLevel'}->{'max_weight'}->[0]->{'ToLevel'}->{'Fixed'}) || # Commented as part of the bug fix of Aug 24 2017
              ($T->{'FromLevel'}->{'max_weight'}->[1] < $max_weight)
            )
          ) {
            $T->{'FromLevel'}->{'Subst'} = [$L,$weight];
          }
        }
      }
    }
  }
  return $num_subst;
} ##FindSubstitutions($$$)

############################################################################
sub GetSubstList($) {  #11/24/2009 6:35PM
############################################################################
  # Recursively build the complete list of substitution levels of the given level.
  my $L = shift;
  my @list = ();
  my $L_subst = defined($L->{'Subst'}) ? $L->{'Subst'}->[0] : undef;
  if ( defined($L_subst) ) {
    push(@list,$L_subst);
    my @list1 = GetSubstList($L_subst);
    push(@list,@list1);
  }
  return @list;
} ##GetSubstList($)

############################################################################
sub GetIndex($$) {   #11/24/2009 4:01PM
############################################################################
  # Return the index of the given level in the current group (with number $GNo).
  my ($L,$GNo) = @_;
  my $k = ($L->{'CurGNo'} == $GNo) ? $L->{'Index'} : -1;
  return $k;
} ##GetIndex($$)

sub BuildMatrixAndVector($$$;$) {
  # Build the elements of the matrix $A and the right-hand parts @E of the
  # least-squares system of equations A*X = E.
  # If $k > 0, add random noise to the values of $A (if $k is odd)
  # or $E (if $k is even).

  my ($G,$GNo,$Pass,$k) = @_;

  # Initialize the matrix $A_terms and vector $E_terms used for
  # bookkeeping of the terms in the sums of which the elements of A and E
  # are constructed.
  my $A_terms = [];
  my $E_terms = [];
  for (my $i = 0; $i < $NA; $i++) {
    $A_terms->[$i] = [];
    $E_terms->[$i] = [];
    for (my $j = 0; $j < $NA; $j++) {
      $A_terms->[$i]->[$j] = [];
    }
  }
           # Main cycle over all transitions
  for (my $i = 1; $i <= $#GroupLines + 1; $i++) {
    my $PT = $GroupLines[$i-1];
    next if ($PT->{'Flags'} & ($FMP | $FExcluded));
    my $PLfrom = $PT->{'FromLevel'};
    my $PLto = $PT->{'ToLevel'};
    if ( ($PLto->{'Name'} eq '110 4') && ($PLfrom->{'Name'} eq '0 0') ||
       ($PLfrom->{'Name'} eq '110 4') && ($PLto->{'Name'} eq '0 0')
    ) {
      $PLto = $PLto;
    }
    my $k1 = $PLto->{'Index'};
    my $k2 = $PLfrom->{'Index'};
    next if (($PLto->{'CurGNo'} != $GNo ) || ($PLfrom->{'CurGNo'} != $GNo));

    my $FF = abs($Sigma[$i]);   # Weight of transition E(k1) -> E(k2) (k1 is the lower level, k2 is the upper level)
    my $Eobs = $PT->{'Eobs'};   # Observed transition wavenumber in cm-1
    my $unc = 1/sqrt($FF);
    if ( $k ) {
      my $delta;
      $trial_size = 1;
      $delta = (2*int(rand(2)) - 1)*$trial_size;
      if ( $k == 2*int($k/2) ) {
        $Eobs += 1.1*$eps*$Eobs*$delta;
      } else {
        $unc += 1.1*$delta*$eps*$unc;
        $FF = 1/($unc*$unc);
      }
    }

    my @substListFrom = ($PLfrom,GetSubstList($PLfrom));
    my @substListTo = ($PLto,GetSubstList($PLto));

    for ( my $j1 = 0; $j1 <= $#substListFrom; $j1++ ) {
      next unless exists $substListFrom[$j1];
      my $L = $substListFrom[$j1];
      my $j2 = IndexOf($L,\@substListTo);
      if ( $j2 >= 0 ) {
        # Delete this level from both lists @substListFrom and @substListTo
        delete($substListFrom[$j1]);
        delete($substListTo[$j2]);
      }
    }
    for ( my $j1 = 0; $j1 <= $#substListTo; $j1++ ) {
      next unless exists $substListTo[$j1];
      my $L = $substListTo[$j1];
      my $j2 = IndexOf($L,\@substListFrom);
      if ( $j2 >= 0 ) {
        # Delete this level from both lists @substListFrom and @substListTo
        delete($substListTo[$j1]);
        delete($substListFrom[$j2]);
      }
    }

    my $hasK21 = 0;
    for (my $j2 = 0; $j2 <= $#substListFrom; $j2++) {
      next unless exists($substListFrom[$j2]);
      my $L2 = $substListFrom[$j2];
      my $k21 = $L2->{'Index'};
      #   A Kramida Aug 25 2017 bug fix for cases when "FromLevel" level of $PT
      #   is the ground state ($Eobs is negative)
      #if (($k21 < $NA) && ($L2->{'CurGNo'} == $GNo)) {
      if ( (($k21 < $NA) || ($j2 == 0)) && ($L2->{'CurGNo'} == $GNo) ) {
      #   End of A Kramida Aug 25 2017 bug fix
        # $j2 and $k21 are indexes of $PLfrom or its substituted levels
        push(@{$A_terms->[$k21]->[$k21]},[$k21,$k21,1,$FF]);
        push(@{$E_terms->[$k21]},[$k21,1,$Eobs*$FF]);
        $hasK21++;
        for (my $j22 = $j2+1; $j22 <= $#substListFrom; $j22++) {
          # $j22 and $k22 are indexes of substituted levels of $PLfrom
          next unless exists($substListFrom[$j22]);
          my $L22 = $substListFrom[$j22];
          my $k22 = $L22->{'Index'};
          if ( ($k22 < $NA) && ($L22->{'CurGNo'} == $GNo) ) {
            push(@{$A_terms->[$k21]->[$k22]},[$k21,$k22,1,$FF]);
            push(@{$A_terms->[$k22]->[$k21]},[$k22,$k21,1,$FF]);
          }
        }
      }
      for (my $j1 = 0; $j1 <= $#substListTo; $j1++) {
        next unless exists($substListTo[$j1]);
        my $L1 = $substListTo[$j1];
        my $k11 = $L1->{'Index'};
        if ( ($k11 < $NA) && ($L1->{'CurGNo'} == $GNo) ) {
          # $j1 and $k11 are indexes of $PLto or its substituted levels
          push(@{$A_terms->[$k11]->[$k21]}, [$k11,$k21,-1,$FF]);
          push(@{$A_terms->[$k21]->[$k11]}, [$k21,$k11,-1,$FF]);
          if ($hasK21 <= 1) {
            push(@{$A_terms->[$k11]->[$k11]}, [$k11,$k11,1,$FF]);
            push(@{$E_terms->[$k11]},[$k11,-1,$Eobs*$FF]);
            for (my $j12 = $j1+1; $j12 <= $#substListTo; $j12++) {
              # $j12 and $k12 are indexes of substituted levels of $PLto
              next unless exists($substListTo[$j12]);
              my $L12 = $substListTo[$j12];
              my $k12 = $L12->{'Index'};
              if ( ($k12 < $NA) && ($L12->{'CurGNo'} == $GNo) ) {
                push(@{$A_terms->[$k11]->[$k12]},[$k11,$k12,1,$FF]);
                push(@{$A_terms->[$k12]->[$k11]},[$k12,$k11,1,$FF]);
              }
            }
          }
        }
      }
    }
    unless ( $hasK21 ) {
      for (my $j1 = 0; $j1 <= $#substListTo; $j1++) {
        next unless exists($substListTo[$j1]);
        my $L1 = $substListTo[$j1];
        my $k11 = $L1->{'Index'};
        if (($k11 < $NA) && ($L1->{'CurGNo'} == $GNo) ) {
          # $j1 and $k11 are indexes of $PLto or its substituted levels
          push(@{$A_terms->[$k11]->[$k11]}, [$k11,$k11,1,$FF]);
          push(@{$E_terms->[$k11]},[$k11,-1,$Eobs*$FF]);
          for (my $j12 = $j1+1; $j12 <= $#substListTo; $j12++) {
            next unless exists($substListTo[$j12]);
            my $L12 = $substListTo[$j12];
            my $k12 = $L12->{'Index'};
            if ( ($k12 < $NA) && ($L12->{'CurGNo'} == $GNo) ) {
              # $j12 and $k12 are indexes of substituted levels of $PLto
              push(@{$A_terms->[$k11]->[$k12]},[$k11,$k12,1,$FF]);
              push(@{$A_terms->[$k12]->[$k11]},[$k12,$k11,1,$FF]);
            }
          }
        }
      }
    }
  }

    # Use the bookkeeping held in $A_terms to calculate the A matrix
  for (my $i = 0; $i < $NA; $i++) {
    for (my $j = 0; $j < $NA; $j++) {
      my $Aij = [];
      my $N_terms = $#{$A_terms->[$i]->[$j]} + 1;
      for ( my $m = 0; $m < $N_terms; $m++ ) {
        my ($k1,$k2,$f,$FF) = @{$A_terms->[$i]->[$j]->[$m]};
        next unless $f;
        my $cancel = 0;
        for ( my $n = $m + 1; $n < $N_terms; $n++ ) {
          my ($k3,$k4,$f1,$FF1) = @{$A_terms->[$i]->[$j]->[$n]};
          if ( #($k3 == $k1) && ($k4 == $k2) && ($f1 == -$f) ||
            ($f1 == -$f) && (abs($FF - $FF1) < 0.5*$eps*abs($FF + $FF1) )
          ) {
            $cancel = 1;
            $A_terms->[$i]->[$j]->[$m]->[2] = 0;
            $A_terms->[$i]->[$j]->[$n]->[2] = 0;
            last;
          }
        }
        unless ($cancel) {
          push(@{$Aij},$FF*$f);
        }
      }
      $A->[$i]->[$j] = Kohan_sum($Aij);
    }
    # Use the bookkeeping held in $E_terms to calculate the E vector
    my $N_terms = $#{$E_terms->[$i]} + 1;
    my $Ei = [];
    for ( my $m = 0; $m < $N_terms; $m++ ) {
      my ($k,$f,$FF) = @{$E_terms->[$i]->[$m]};
      next unless $f;
      my $cancel = 0;
      for ( my $n = $m + 1; $n < $N_terms; $n++ ) {
        my ($k1,$f1,$FF1) = @{$E_terms->[$i]->[$n]};
        if ( #($k1 == $k) && ($f1 == -$f) ||
          ($f1 == -$f) && (abs($FF - $FF1) < 0.5*$eps*abs($FF + $FF1))
        ) {
          $cancel = 1;
          $E_terms->[$i]->[$m]->[1] = 0;
          $E_terms->[$i]->[$n]->[1] = 0;
          last;
        }
      }
      unless ($cancel) {
        push(@{$Ei}, $FF*$f);
      }
    }
    $E[$i] = Kohan_sum($Ei);
  }
}

my $scale_row = [];
my $scale_col = [];

sub FormMatrix($$$) {
  # Parameters: Level group $G, its index $GNo, and the stage of optimization $Pass.
  # Actions:
  # - Add all fixed levels to the current group.
  # - Reorder $G so that the ground level is the last in the list.
  # - Initialize and populate the array of transition weights $Sigma.
  # - Exclude transitions that cannot be processed with this group of levels.
  # - For each level, find variable substitutions if necessary (by calling FindSubstitutions()).
  # - For each level, count the numbers of connecting transitions of various categories.
  # - Initialize and populate the matrix A and the vector of right-hand parts of the
  #   system of equations A*X = E (by calling BuildMatrixAndVector())

  my ($G,$GNo,$Pass) = @_;

  my ($k,$NL);

  $numFixLevs = 0;
  # Bug fix A. Kramida Aug 24 2017
  my %hash = ();

  $hash{$Ground} = $Ground;
  foreach my $L (@{$G}) {
    next if defined $hash{$L};
    foreach my $PT (@VirtLines) {
      next if ( $PT->{'FromLevel'} eq $PT->{'ToLevel'} );
      my $L1 = $PT->{'FromLevel'};
      my $L2 = $PT->{'ToLevel'};
      next if ( ($L1 ne $L) && ($L2 ne $L) );
      $hash{$L1} = $L1;
      $hash{$L2} = $L2;
    }
  }
  my @fixedThisGroup = values %hash;

  #foreach my $FL (@FixLevs) {
  foreach my $FL (@fixedThisGroup) {
  # End bug fix A. Kramida Aug 24 2017
    # Insert all fixed levels into the group
    InsertFixedLev($FL,$G);
    $numFixLevs++ if ($FL->{'Fixed'} > 0);
  }

  # Determine the dimension of the matrix, $NA
  $NA = $#{$G}; # One less than level count in the group

  my $i = IndexOf($Ground,$G);
  if ($i < 0) {
    Error('Ground level not found in group $GNo.');
  }
        # Make the ground level the last in the group
  my $PL = $G->[$NA];
  $G->[$NA] = $Ground;
  $G->[$i] = $PL;

  return if ($NA==0);

  $num_subst = FindSubstitutions($G,$Pass,$num_subst);

  if ( ($Pass == 2) && ($numFixLevs == 1) && !$num_subst) {
    # If there is only one truly fixed level, i.e. the ground level,
    # and there were no variable substitutions affecting the matrix,
    # there is no need to re-solve the matrix inversion problem,
    # since the matrix is going to be the same as in Pass 1.
    return;
  }

  $A = [];
  for (my $i = 0; $i < $NA; $i++) {
    $E[$i] = 0;
    $A->[$i] = [];
    for ( my $j = 0; $j < $NA; $j++ ) {
      $A->[$i]->[$j] = 0;
    }
    $Sigma[$i] = 0;
  }

  $k = 0;
  if ($Pass==1) {
    @GroupLines = ();
    $Fa = [];
    $fb = [];
    $r = [];
    for (my $i = 0; $i < $NA; $i++) {
      $Fa->[$i] = [];
    }
    my $i1 = 0;
    foreach my $Level (@{$G}) {
      $i1++;
      $Level->{'NV'} = 0;
      $Level->{'CurGNo'} = $GNo;
      $Level->{'Index'} = $i1-1;
    }
  }

  $NL = ($Pass==2) ? $#GroupLines + 1 : $#Lines + 1;   # Number of lines
  for (my $i = 1; $i <= $NL; $i++) {
    my $PT = ($Pass==1) ? $Lines[$i-1] : $GroupLines[$i-1];
    $PT->{'InGroup'} = ($Pass==1) ? 0 : $i;
    next if ($PT->{'Flags'} & $FExcluded);

    my $j1 = $PT->{'FromLevel'}->{'Index'};
    my $j2 = $PT->{'ToLevel'}->{'Index'};

    if (($PT->{'FromLevel'}->{'CurGNo'} != $GNo) || ($PT->{'ToLevel'}->{'CurGNo'} != $GNo)) {
      if ( (($PT->{'FromLevel'}->{'CurGNo'} == $GNo) || ($PT->{'ToLevel'}->{'CurGNo'} == $GNo)) && ($PT->{'Flags'} & $FMP) ) {
          # Predicted or masked intersystem line
        if (!$PT->{'ToLevel'}->{'Fixed'} && !$PT->{'FromLevel'}->{'Fixed'}) {
          print 'Set "Divide level groups=N" ' ,
            "in par. file to predict line No. $i (intersystem).\n";

          $PT->{'Flags'} |= $FExcluded;
          $PT->{'InGroup'} = 0;
        }
      }
      next;
    }
    $k++;
    if ($Pass==1) {
      # Count the numbers of connecting transitions of various categories
      $PT->{'InGroup'} = $k;
      push(@GroupLines,$PT);

      if ($PT->{'Flags'} & $FMP) {
        $PT->{'FromLevel'}->{'NM'}++;
        $PT->{'ToLevel'}->{'NM'}++;
        next;
      }
      if ($PT->{'Flags'} & $FVirtual) {
        if ($PT->{'WnUnc'} < $too_small) {
          Error('Zero uncertainty is given for more than one fixed level.');
        }
        $PT->{'FromLevel'}->{'NV'}++;
        $PT->{'ToLevel'}->{'NV'}++;
      } else {
        Error('Error: zero line uncertainty.') if ($PT->{'WnUnc'} < $too_small);
        $PT->{'FromLevel'}->{'NT1'}++;
        $PT->{'ToLevel'}->{'NT1'}++;
      }
    }

    $Sigma[$k] = GetLineWeight($PT,$Pass);
  }
  $NLines = $k;

  &BuildMatrixAndVector($G,$GNo,$Pass);

}  # FormMatrix

sub FillD($$$$$$) {
  # Determine the various uncertainty values for level $Lev
  # from the elements of inverted matrix B = A^-1.
  # (Called by FindLevels())

  my ($Lev,$i,$NA,$G1,$PLg,$GroupNo) = @_;
  my ($EL,$FF);

  $i++;
  return $i if ($i>$NA);
  my $D = $B->[$i-1]->[$i-1];

  $D = ($D > 0) ? sqrt($D) : 0;

  my $IsFixed = $Lev->{'Fixed'};
  if ( !$IsFixed || ($D > $Lev->{'D1'}) ) {
    $Lev->{'D1'} = $D;
  }

  my $k = $PLg->{'Index'};
  my $PLgIsGround = ($PLg eq $Ground);
  if ( ($PLg->{'CurGNo'} == $GroupNo) || $PLgIsGround ) {
    if ($PLgIsGround) {
      $D = 0;
    } else {
      $D = abs($B->[$i-1]->[$i-1] + $B->[$k]->[$k] - 2*$B->[$i-1]->[$k]);
    }
    $D = sqrt($D);
    if (!$IsFixed || (($D > $Lev->{'D2'}) && ($Lev->{'GNo'} == $GroupNo))) {
      $Lev->{'D2'} = $D;
    }
  }
    # Find uncertainties due to possible systematic shifts of
    #           correlated line groups
  if ( !$IsFixed ) {
    for (my $jj = 0; $jj < $MaxCGNo; $jj++) {
      my $PCG = $CGroups[$jj];
      for (my $k = 1; $k <= $#{$PCG} + 1; $k++) {
        my $PT1 = $PCG->[$k-1];

        my $m = $PT1->{'InGroup'};
        next if ($m <= 0);

        my $k1 = $PT1->{'FromLevel'}->{'Index'};
        next if ($PT1->{'FromLevel'}->{'CurGNo'} != $GroupNo);
        my $To = $PT1->{'ToLevel'};
        my $k2 = $To->{'Index'};
        next if ($To->{'CurGNo'} != $GroupNo);

        $D  = -$B->[$k1]->[$i-1];
        if ($To ne $Ground) {
          $D += $B->[$k2]->[$i-1];
        }
        $Lev->{'C'}->[$jj] += $D*abs($Sigma[$m])*$CUnc->[$jj];
      }
    }
  }
  return $i;
}

############################################################################
sub SubstBack($) {  #11/25/2009 2:56PM
############################################################################
  # Determine the level energies from the substituted variables and
  # their values found in the matrix inversion procedure with iterative
  # correction.
  my ($G) = @_;
  foreach my $L (@{$G}) {
    next if defined($L->{'Esubst'});
    my @substLevs = GetSubstList($L);
    my $num_s = $#substLevs;
    if ( $num_s < 0) {
      $L->{'Esubst'} = $L->{'E'};
    } else {
      my $E = 0;
      for (my $i = $num_s; $i >= 0; $i--) {
        my $L1 = $substLevs[$i];
        if ( $i == $num_s ) {
          $E = $L1->{'E'};
        } else {
          $E += $L1->{'E'};
        }
        $L1->{'Esubst'} = $E;
      }
      $L->{'Esubst'} = $L->{'E'} + $E;
    }
  }
  foreach my $L (@{$G}) {
    if ( defined($L->{'Esubst'}) ) {
      $L->{'E'} = $L->{'Esubst'};
      $L->{'Esubst'} = undef;
    }
  }
} ##SubstBack($)

############################################################################
sub SubstBack1($$$) {  #11/25/2009 2:56PM
############################################################################
  # For random error trials,
  # Determine the level energies from the substituted variables and
  # their values found in the matrix inversion procedure with iterative
  # correction.
  my ($G,$n,$x) = @_;
  my $subst = [];
  for (my $j = 0; $j < $n; $j++) {
    my $L = $G->[$j];
    next if defined($subst->[$j]);
    my @substLevs = GetSubstList($L);
    my $num_s = $#substLevs;
    if ( $num_s < 0) {
      $subst->[$j] = $x->[$j];
    } else {
      my $E = 0;
      for (my $i = $num_s; $i >= 0; $i--) {
        my $L1 = $substLevs[$i];
        my $k = $L1->{'Index'};
        if ( $i == $num_s ) {
          $E = $x->[$k];
        } else {
          $E += $x->[$k];
        }
        $subst->[$k] = $E;
      }
      $subst->[$j] = $x->[$j] + $E;
    }
  }
  return $subst;
} ##SubstBack1($$$)

############################################################################
sub scale_rows($) {  #04/23/2010 11:22AM
############################################################################
  my $A = shift;
  my $A_scaled = [];
  $scale_row = [];
  if ( $use_Cholesky ) {
    for ( my $m = 0; $m < $NA; $m++ ) {
      $scale_row->[$m] = 1;  # No row scaling for Cholesky
      for ( my $n = 0; $n < $NA; $n++ ) {
        $A_scaled->[$m]->[$n] = $A->[$m]->[$n];
      }
    }
  } else {
    # Scale rows to minimize computational errors in matrix inversion
    for ( my $m = 0; $m < $NA; $m++ ) {
      my $row_sum = 0;
      for ( my $n = 0; $n < $NA; $n++ ) {
        $row_sum += abs($A->[$m]->[$n]);
      }
      $scale_row->[$m] = ($row_sum > $too_small) ? 1/$row_sum : 1;
      for ( my $n = 0; $n < $NA; $n++ ) {
        $A_scaled->[$m]->[$n] = $A->[$m]->[$n]*$scale_row->[$m];
      }
    }
  }
  return $A_scaled;
} ##scale_rows($)

############################################################################
sub scale_cols($) {  #05/14/2010 8:25AM
############################################################################
  # Start with the row-scaled matrix.
  my $A_scaled = shift;
  $scale_col = [];
      # Scale columns to minimize computational errors in matrix inversion.
  for ( my $n = 0; $n < $NA; $n++ ) {
    my $col_sum = 0;
    for ( my $m = 0; $m < $NA; $m++ ) {
      $col_sum += abs($A_scaled->[$m]->[$n]);
    }
    $scale_col->[$n] = $col_sum > $too_small ? 1/$col_sum : 1;
    for ( my $m = 0; $m < $NA; $m++ ) {
      $A_scaled->[$m]->[$n] *= $scale_col->[$n];
    }
  }
  return $A_scaled;
} ##scale_cols($)

sub FindLevels($$$) {
  # Parameters: Levels group $Group, its index $GroupNo, stage of optimization $Pass.
  # Actions:
  # - Determine if matrix inversion is needed (always yes for Pass==1).
  # - If yes, then populate the matrix-scaling vector $scale_row,
  #   scale the A matrix to A_scaled, invert it to B,
  #   scale the right-hand part vector @E to $b,
  #   find the solution to A_scaled * x = b,
  #   if Pass == 1 then
  #      1) make an iterative improvement of the solution x,
  #      2) find the numerical errors in x and scale them to errors in the level energies,
  #      3) find the matrix condition number,
  #      4) find the lowest level in the group $PLg.
  #   Scale back the matrix B.
  # - If Pass == 2 then
  #      1) determine the various uncertainty values for each level by calling FillD(),
  #      2) sum up deviations of observed wavenumbers from their Ritz values to determine
  #         the Radziemski uncertainty D for each level (designated D1 in the levels output file).

  my ($Group,$GroupNo,$Pass) = @_;
  my $i = 0;

  my ($rms,$big_dev);
  if ( $Pass == 1 ) {
    my $xs = [];
    my $xd = [];
    my $yd = [];
    for ( my $k = 0; $k <= 2*$statistics_size; $k++ ) {
      my $x = [];
      &BuildMatrixAndVector($Group,$GroupNo,$Pass,$k) if $k;

      my $A_scaled = &scale_rows($A);

      # Scale columns to minimize computational errors in matrix inversion.
      #  This scaling is turned off in the present version.
      #$A_scaled = &scale_cols($A_scaled);

         #  Invert the scaled matrix
      my ($L,$piv) = $use_Cholesky ? CholeskyDecomposition($A_scaled) : LUDecomposition1($A_scaled);

      if ( ref($L) ne 'ARRAY' ) {
        if ( $L < 0 ) {
          print "\n   Matrix is asymmetric!\n";
        } else {
          print "\n   Degenerate matrix for energies. Try to fix the level:\n",
            $Group->[$L]->{'Name'},"\n";
        }
        &PrintMatrix($Group,$NA);
        exit;
      }

      #$cond = &FindConditionNumber($A_scaled,$B);
      #$rel_cond = &FindRelCondNumber($A_scaled,$B);
      #$normwise_error_bound = $eps*$cond;
      #$normwise_error_bound /= (1-$eps*$cond/2) if ($eps*$cond/2 < 1);

      my $b = [];
      for ( my $m = 0; $m < $NA; $m++ ) {
        $b->[$m] = [$scale_row->[$m]*$E[$m]];
      }

         # Find the scaled solution $x
      $x = $use_Cholesky ? SolveCholesky($L,$b) : LUSolve($L,$b,$piv);

      for ( my $m = 0; $m < $NA; $m++ ) {
        $x->[$m] = $x->[$m]->[0];
      }
         # Find residuals
      $r = [];
      my $max_r = 0;
      my $max_Ax = 0;
      my $min_Ax = 1e100;
      for ( my $m = 0; $m < $NA; $m++ ) {
        my $Ax = 0;
        my $rm = [-$b->[$m]->[0]];
        for ( my $n = 0; $n < $NA; $n++ ) {
          push(@{$rm}, $A_scaled->[$m]->[$n]*$x->[$n]);
          $Ax += abs($A_scaled->[$m]->[$n]*$x->[$n]);
        }
        $r->[$m] = Kohan_sum($rm);
        $max_r = abs($r->[$m]) if $max_r < abs($r->[$m]);
        $max_Ax = $Ax if $max_Ax < $Ax;
        $min_Ax = $Ax if ($Ax && ($min_Ax > $Ax));
      }

      #$iter_crit = $min_Ax ? $cond * $max_Ax/$min_Ax * $eps : 0;

      my $iterate = $use_Cholesky ? 0 : 1; #($iter_crit < 1) ? 1 : 0;

      if ( $iterate ) {
        my $d = [];
        for ( my $m = 0; $m < $NA; $m++ ) {
          $d->[$m] = [$r->[$m]];
        }
          # Find an iterative refinement
        $d = $use_Cholesky ? SolveCholesky($L,$d) : LUSolve($L,$d,$piv);
        my $max_d = 0;
        for ( my $m = 0; $m < $NA; $m++ ) {
          $d->[$m] = $d->[$m]->[0];
          $max_d = $d->[$m] if ($max_d < $d->[$m]);
        }
        # Apply the iterative correction
        for ( my $m = 0; $m < $NA; $m++ ) {
          $x->[$m] -= $d->[$m];
        }
        # Find new residuals
        $r = [];
        $max_r = 0;
        for ( my $m = 0; $m < $NA; $m++ ) {
          my $rm = [-$b->[$m]->[0]];
          for ( my $n = 0; $n < $NA; $n++ ) {
            push(@{$rm}, $A_scaled->[$m]->[$n]*$x->[$n]);
          }
          $r->[$m] = Kohan_sum($rm);
          $max_r = abs($r->[$m]) if $max_r < abs($r->[$m]);
        }
      }

      if ( $k == 0 ) {
        for ( my $m = 0; $m < $NA; $m++ ) {
          $xs->[$m] = $x->[$m];
          my $PL = $Group->[$m];
          $PL->{'E'} = $x->[$m];
        }
        SubstBack($Group) if $num_subst;

          # Find the lowest level in the group, $PLg
        $PLg = undef;
        foreach my $PL (@{$Group}) {
          if ($PL->{'Fixed'} && ($PL->{'GNo'} == $GroupNo)
           && (!defined($PLg) || ($PL->{'E'} < $PLg->{'E'}))
          ) {
           $PLg = $PL;
          }
        }
        $PLg = $Ground unless (defined($PLg));
      } else {
        my $xx = SubstBack1($Group,$NA,$x);
        for ( my $m = 0; $m < $NA; $m++ ) {
          my $d = $xx->[$m] - $Group->[$m]->{'E'};
          $d *= $d;
            # Accumulate the sum of squares of level deviations
          if ( $k == 2*int($k/2) ) {
            $yd->[$m] += $d;  # Accumulated sum of squares of errors due to small changes in Eobs of transitions
          } else {
            $xd->[$m] += $d;  # Accumulated sum of squares of errors due to small changes in measurement uncertainties of transitions
          }
        }
      }
    }

    if ( $statistics_size ) {
      my $max_d = 0;
      my $max_d1 = 0;
      for ( my $m = 0; $m < $NA; $m++ ) {
        my $PL = $Group->[$m];
        $PL->{'d'} = 0.21*sqrt(($xd->[$m]+$yd->[$m])/$statistics_size/$trial_size);
        $max_d = $PL->{'d'} if $max_d < $PL->{'d'};
      }
      $max_d = $max_d;
    }
  } else {
      # The 3rd parameter of InvertMatrix below is 0 for LU decomposition and 1 for Cholesky.
      # Testing showed that inversion with LU is faster than Cholesky in this Perl version.
      # Since the numerical stability is not critical at this stage of the computation,
      # I use LU decomp and LUSolve to invert the matrix for finding the level and transition uncertainties.
    ($B,$rms,$big_dev) = InvertMatrix($A,$InvChk,0);
    if ( $B < 0 ) {
      print "\n   Matrix is degenerate!\n";
      &PrintMatrix($Group,$NA);
      exit;
    }

    $i = 0;
    foreach my $PL (@{$Group}) {
      $i = FillD($PL,$i,$NA,$Group,$PLg,$GroupNo);
    }

    # Sum up line deviations for each level

    for (my $i = 1; $i <= $NLines; $i++) {
      my $PT = $GroupLines[$i-1];
      my $From = $PT->{'FromLevel'};
      next if ( $From->{'GNo'} != $GroupNo );
      my $To = $PT->{'ToLevel'};
      next if ( $To->{'GNo'}   != $GroupNo );

      my $dE = $PT->{'WnUnc'};
      if ( ($dE > $too_small) and !($PT->{'Flags'} & $FMP) ) {
        my $dE1 = trans_Ecalc($PT) - $PT->{'Eobs'};
              # Add to each level the weighted line deviation
              # and the weighted experimental line uncertainty
        my $FF = 1/sqr($dE);

        if ( $FF > $too_small ) {
          $From->{'NT'} += $FF;  # Sum of all weights
          $To->{'NT'}   += $FF;
          my $dd = $FF*(1+sqr($dE1/$dE));
          $From->{'D'} += $dd;
          $To->{'D'}   += $dd;
        }
      }
    }
  }
}

sub FindUncert($$) {
          # Calculate uncertainties of predicted wavenumbers
  my ($G1,$GNo) = @_;

  for (my $i = 1; $i <= $#GroupLines + 1; $i++) {
    my $PT = $GroupLines[$i-1];
    my $PL2 = $PT->{'FromLevel'};
    my $PL  = $PT->{'ToLevel'};

    my $i1 = $PL2->{'Index'};
    next if ($PL2->{'CurGNo'} != $GNo);
    my $i2 = $PL->{'Index'};
    next if ($PL->{'CurGNo'} != $GNo);

    my $D = ($PL eq $Ground) ? abs($B->[$i1]->[$i1]) : ($PL2 eq $Ground) ?
      abs($B->[$i2]->[$i2]) :
      abs($B->[$i1]->[$i1] + $B->[$i2]->[$i2] - 2*$B->[$i1]->[$i2]);

    $PT->{'WnUncCalc'} = sqrt($D);

      # Find uncertainties due to possible systematic shifts of
      #           correlated line groups
    next if (($PT->{'Flags'} & $FVirtual) && !($PT->{'Flags'} & $FPredict));

    my $Dc = 0;
    my $PLnotGroundNotPLg = (($PL ne $Ground) && ($PL ne $PLg));
    for (my $j = 0; $j < $MaxCGNo; $j++) {
      my $PCG = $CGroups[$j];
      for (my $k = 0; $k <= $#{$PCG}; $k++) {
        my $PT1 = $PCG->[$k];

        my $m = $PT1->{'InGroup'};
        next if ($m <= 0);

        my $From1 = $PT1->{'FromLevel'};
        my $k1 = $From1->{'Index'};
        next if ($From1->{'CurGNo'} != $GNo);
        my $To1 = $PT1->{'ToLevel'};
        my $k2 = $To1->{'Index'};
        next if ($To1->{'CurGNo'} != $GNo);

        my $FF = -$B->[$k1]->[$i1];
        if ($PLnotGroundNotPLg) {
          $FF += $B->[$k1]->[$i2];
        }
        if ($PT1->{'ToLevel'} ne $Ground) {
          $FF += $B->[$k2]->[$i1];
          if ($PLnotGroundNotPLg) {
            $FF -= $B->[$k2]->[$i2];
          }
        }
        $PT->{'C'}->[$j] += $FF*abs($Sigma[$m]) * $CUnc->[$j];
      }
      $Dc += sqr($PT->{'C'}->[$j]);
    }

    if (trans_SingleLineLevelFlag($PT) ne '*') {
      $Dc = sqrt($Dc + sqr($PT->{'WnUncCalc'}));
      if ( ($Dc <= $PT->{'WnUnc'}) || ($PT->{'Flags'} & $FMP) ) {
        # Increase the calculated line uncertainty due to line correlations
        # but assume that it cannot exceed the measured line uncertainty.
        # Yet, if the calculated (random) dispersion already exceeds the
        #  declared measured uncertainty, leave it so
        $PT->{'WnUncCalc'} = $Dc;
      } elsif ($PT->{'WnUncCalc'} < $PT->{'WnUnc'}) {
        $PT->{'WnUncCalc'} = $PT->{'WnUnc'};
      }

      if ( $PT->{'Flags'} & ($FVirtual | $FPredict) ) {
        # Determine which of the two levels is the base one
        my ($BL,$L) = (undef,undef);
        my $k1 = IndexOf($PL,\@UncBases);

        if ( $k1 < 0 ) {
          $k1 = IndexOf($PL2,\@UncBases);
          if ($k1 >= 0) {
            $BL = $PL2;
            $L = $PL;
          }
        } else {
          $BL = $PL;
          $L = $PL2;
        }

        if ( defined $BL ) {
          # Set the value of uncertainty relative to the base level
          $L->{'UncRelBases'}->[$k1] = $PT->{'WnUncCalc'};
          $BL->{'UncRelBases'}->[$k1] = '_';
        }
      }
    }
  }
    # Find max fraction of numerical errors in the level energies
  $max_num_error_frac = 0;
  foreach my $L (@{$G1}) {
    my $num_error = $L->{'d'};
    my $min_uncert = $L->{'D'};
    # If the level is determined by only one line or is fixed,
    # ignore the small $L->{'D'} and replace it with a larger $L->{'D1'} or $L->{'D2'},
    # since adjustment of the single-line levels guarantees the exactness
    # of the corresponding Ritz wavenumber.
    $min_uncert = $L->{'D1'} if ($min_uncert > $L->{'D1'} || $L->{'NT2'} <= 1);
    $min_uncert = $L->{'D2'} if (($L->{'D2'} > 0) && $min_uncert > $L->{'D2'});
    foreach my $d (@{$L->{'UncRelBases'}}) {
      $min_uncert = $d if ($d < $min_uncert) && ($d > 0);
    }
    my $frac = $min_uncert ? $num_error/$min_uncert : 0;
    $max_num_error_frac = $frac if $max_num_error_frac < $frac;
  }
} # FindUncert

sub FindWidths() {
  # Calculate widths of columns in the output tables.

  for (my $i = 1; $i <= $#Levels + 1; $i++) {
     # Find widths of columns in the levels output file
    my $PL = $Levels[$i-1];

    foreach my $s (qw{E D D1 D2 d}) {
      my $decPl = num_dec_places($PL->{"${s}r"});
      $Widths{"pl$s"} = $decPl if ($decPl > $Widths{"pl$s"});
      my $w = length($PL->{"${s}r"});
      $w -= ($decPl + 1) if ($decPl > 0);
      $Widths{"width$s"} = $w if ($w > $Widths{"width$s"});
    }
    my $LName_length = length($PL->{'Name'});
    $Widths{'LName'} = $LName_length if $LName_length > $Widths{'LName'};

    for ( my $i = 0; $i <= $#UncBases; $i++ ) {
      unless (defined $Widths{'Bases'}) {
        $Widths{'Bases'} = [];
        $Widths{'BasePlaces'} = [];
      }
      my $decPl = num_dec_places($PL->{'UncRelBases'}->[$i]);
      $Widths{'BasePlaces'}->[$i] = $decPl
        if ($decPl > $Widths{'BasePlaces'}->[$i]);
      my $w = length($PL->{'UncRelBases'}->[$i]);
      $w -= ($decPl + 1) if ($decPl > 0);
      $Widths{'Bases'}->[$i] = $w if ($w > $Widths{'Bases'}->[$i]);
    }
  }
  for ( my $i = 0; $i <= $#UncBases; $i++ ) {
    $Widths{"Bases"}->[$i] += ($Widths{"BasePlaces"}->[$i] + ($Widths{"BasePlaces"}->[$i] > 0 ? 1 : 0) + 1);
  }

  for (my $i = 1; $i <= $#Lines + 1; $i++) {
     # Find widths of columns in the lines output file
    my $PT = $Lines[$i-1];
    next if ( ($PT->{'Flags'} & $FVirtual) && !$WriteVirt);
    next if ($PT->{'Flags'} & $FExcluded);

    foreach (qw{_Lobs _Eobs _LuncObs _wnObs _WnUncObs _LairCalc _LvacCalc
      _LuncCalc _WnCalc _WnUncCalc _dW _dE _flags}
    ) {
      my $X = $PT->{$_};
      my $w = length($X);
      my $decimalPlaces = num_dec_places($X);
      $w -= $decimalPlaces;

      $Widths{"width$_"} = $w if $w > $Widths{"width$_"};
      $Widths{"pl$_"} = $decimalPlaces if $decimalPlaces > $Widths{"pl$_"};
    }
    my $I_length = length($PT->{'Iobs'});
    $Widths{'width_Iobs'} = $I_length if $I_length > $Widths{'width_Iobs'};
  };

  foreach (qw{_Lobs _Eobs _LuncObs _WnUncObs _LairCalc _LvacCalc
    _LuncCalc _WnCalc _WnUncCalc _dW _dE E D D1 D2 d}
  ) {
    $Widths{"width$_"} += ($Widths{"pl$_"} + ($Widths{"pl$_"} > 0 ? 1 : 0) + 1);
  }
}

my $SLGroups = [];
my $SLGroupLines = [];
my $SLBaseLevels = [];
my $NeedAdjustment = 1;

############################################################################
sub FillSLGroups {   #10/20/2009 8:45AM
############################################################################
  # Populate groups of single-line levels - this is needed for adjustment of their energies.
  # Set the global variable $NeedAdjustment to 1 if adjustment is needed.
  $SLGroups = [];
  $SLGroupLines = [];
  $SLBaseLevels = [];
  for (my $i=0; $i <= $#SLLines; $i++) {
    my $PT = $SLLines[$i];
    my ($G1,$G2,$PL1,$PL2,$PT1) = (undef,undef,undef,undef,undef);
    my $PL = $PT->{'FromLevel'};

    $PL1 = $PT->{'ToLevel'};
    if ( $PL->{'Adjusted'} ) {
      if ( $PL1->{'Adjusted'} ) {
        next;
      } else {
        $PL1->{'NextToAdjust'} = 1 unless ((($PL1->{'NT1'}) && ($PL1->{'NT2'} > 1)) || $PL1->{'Fixed'});
      }
    } else {
      if ( $PL1->{'Adjusted'} ) {
        $PL->{'NextToAdjust'} = 1 unless ((($PL->{'NT1'}) && ($PL->{'NT2'} > 1)) || $PL->{'Fixed'});
      } else {
        if ( (($PL->{'NT1'}) && ($PL->{'NT2'} > 1)) || $PL->{'Fixed'}) {
          $PL1->{'NextToAdjust'} = 1;
        } elsif ( (($PL1->{'NT1'}) && ($PL1->{'NT2'} > 1)) || $PL1->{'Fixed'}) {
          $PL->{'NextToAdjust'} = 1;
        }
      }
    }
    $PL1 = undef;

    # Try to find lower and upper levels in already existing SL groups
    foreach my $Gx (@{$SLGroups}) {
      foreach my $PLe (@{$Gx}) {
        if (SameLevels($PLe,$PL)) {
          $PL1 = $PLe;
          last;
        }
      }
      if ($PL1) {
        if ($PL1 ne $PL) {
          $PT->{'FromLevel'} = $PL1;
        }
        $G1 = $Gx;
        last;
      }
    }

    $PL = $PT->{'ToLevel'};
    foreach my $Gx (@{$SLGroups}) {
      foreach my $PLe (@{$Gx}) {
        if (SameLevels($PLe,$PL)) {
          $PL2 = $PLe;
          last;
        }
      }
      if ($PL2) {
        if ($PL2 ne $PL) {
          $PT->{'ToLevel'} = $PL2;
        }
        $G2 = $Gx;
        last;
      }
    }

    if ($G1) {
      if ($G2) {
        my $i_found = 0;
        $PT1 = undef;
        if ($G1 ne $G2) {
          $G1 = MergeGroups($G1,$G2,$SLGroups,0);
        }
      } else {
        $PL2 = $PT->{'ToLevel'};
        push(@{$G1}, $PL2) unless ($PL2->{'Adjusted'} || (($PL2->{'NT1'}) && ($PL2->{'NT2'} > 1)) || $PL2->{'Fixed'});
      }
    } else {
      if ($G2) {
        $PL1 = $PT->{'FromLevel'};
        push(@{$G2},$PL1) unless ($PL1->{'Adjusted'} || (($PL1->{'NT1'}) && ($PL1->{'NT2'} > 1)) || $PL1->{'Fixed'});
      } else {
        $PL1 = $PT->{'FromLevel'};
        $PL2 = $PT->{'ToLevel'};
        push(@{$G1}, $PL2) unless ($PL2->{'Adjusted'} || (($PL2->{'NT1'}) && ($PL2->{'NT2'} > 1)) || $PL2->{'Fixed'});
        push(@{$G1}, $PL1) unless ($PL1->{'Adjusted'} || (($PL1->{'NT1'}) && ($PL1->{'NT2'} > 1)) || $PL1->{'Fixed'});
        push(@{$SLGroups},$G1);
      }
    }
  }

  # Populate the SLGroupLines array
  for (my $i = 0; $i <= $#SLLines; $i++) {
    my $PT = $SLLines[$i];
    for ( my $j = 0; $j <= $#{$SLGroups}; $j++ ) {
      my $k = IndexOf($PT->{'ToLevel'},$SLGroups->[$j]);
      if ( $k >= 0 ) {
        $SLGroupLines->[$j] = [] unless defined($SLGroupLines->[$j]);
        push(@{$SLGroupLines->[$j]},$PT);
        last;
      }
      $k = IndexOf($PT->{'FromLevel'},$SLGroups->[$j]);
      if ( $k >= 0 ) {
        $SLGroupLines->[$j] = [] unless defined($SLGroupLines->[$j]);
        push(@{$SLGroupLines->[$j]},$PT);
        last;
      }
    }
  }

  # Find the base level for each SL group
  for (my $i = 0; $i <= $#{$SLGroups}; $i++) {
    my $G = $SLGroups->[$i];
    my $BL = undef;
    foreach my $PL (@{$G}) {
      if ( ((($PL->{'NT1'}) && ($PL->{'NT2'} > 1)) || $PL->{'Fixed'} || $PL->{'NextToAdjust'}) && !$PL->{'Adjusted'} ) {
        $BL = $PL;
        $PL->{'Adjusted'} = 1;
        $NeedAdjustment = 1;
        last;
      }
    }
    $SLBaseLevels->[$i] = $BL if $BL;
  }
} ##FillSLGroups

sub AdjustSingleLineLevels {

  # Populate the @SLLines array
  for (my $i = 0; $i <= $#Lines; $i++) {
    my $PT = $Lines[$i];
    next unless trans_SingleLineLevelFlag($PT) eq '*';
    next if $PT->{'Flags'} & $FVirtual;
    next if abs($PT->{'Eobs'} < $too_small);
    next if $PT->{'Flags'} & ($FMP | $FExcluded);
    $PT->{'FromLevel'}->{'Shift'} = 0;
    $PT->{'FromLevel'}->{'Adjusted'} = 0;
    $PT->{'FromLevel'}->{'NextToAdjust'} = 0;
    $PT->{'ToLevel'}->{'Shift'} = 0;
    $PT->{'ToLevel'}->{'Adjusted'} = 0;
    $PT->{'ToLevel'}->{'NextToAdjust'} = 0;
    push(@SLLines,$PT);
  }

    # Cycle and apply the adjustment by setting the {'Shift'} values
    # of connecting levels for each SL level
    # until no further adjustment is needed.
  while ( $NeedAdjustment ) {
    $NeedAdjustment = 0;
    FillSLGroups();
    if ( $NeedAdjustment ) {
      for (my $i = 0; $i <= $#{$SLGroups}; $i++) {
        my $G = $SLGroups->[$i];
        my $BL = $SLBaseLevels->[$i];
        if ( $BL ) {
          my $E = $BL->{'E'} + $BL->{'Shift'};
          my $E_rounded = $BL->{'Er'};
          my $dE = $E_rounded - $E;
          foreach my $PL (@{$G}) {
            next if ( (($PL->{'NT1'}) && ($PL->{'NT2'} > 1)) || $PL->{'Fixed'} );
            $PL->{'Shift'} += $dE;
            &RoundLevel($PL);
          }
        }
      }
    }
  }
}

############################################################################
sub lev_GetMinUncKey($)  {  #10/20/2009 12:41PM
############################################################################
  # Return the key to the field containing the minimum non-zero uncertainty of a level.
  my $Lev = shift;
  my $min_unc_key = 'D';       # if $Lev->{'D'}==0 then choose $Lev->{'D1'}
  if (($Lev->{$min_unc_key} eq '_') || ($Lev->{'D1'} > 0) && ($Lev->{'D1'} < $Lev->{$min_unc_key})) {
    $min_unc_key = 'D1';
  }
  if (($Lev->{'D2'} > 0) && ($Lev->{'D2'} < $Lev->{$min_unc_key})) {
    $min_unc_key = 'D2';
  }
  return $min_unc_key;
} ##lev_GetMinUncKey($)

############################################################################
sub lev_GetMinUnc($)  {  #10/20/2009 12:41PM
############################################################################
  # Return the minimum non-zero uncertainty value of a level
  # (from the set D, D1, D2).
  my $Lev = shift;
  return $Lev->{lev_GetMinUncKey($Lev)};
} ##lev_GetMinUnc($)

my $NeedTuning = 1;

############################################################################
sub trans_tune_SL_levels($) {    #11/10/2009 3:05PM
############################################################################
  # Tune the rounding precision of a single-line level defined by the
  # given transition, so that the transition's Ritz value of the observed
  # quantity exactly match the observed quantity value.
  # Sets the global variable $NeedTuning to 1 if further tuning is needed.
  my $PT = shift;
  my ($Etr1,$WnUnc) = trans_Ecalc_rounded($PT);
  return 0 if (abs($Etr1) < $too_small);
  my $lambda = trans_ObsQuantityValue($PT);
  my ($lambda1,$u,$lsd) = trans_CalcValueOfObsQuantityRounded($PT);
  $lsd *= 0.5;
  return 0 if (abs($lambda1-$lambda) < $lsd);

    # Add one more significant figure to one of the rounded energy levels.
    # First, try the level having the largest rounding error.
  my $PL1 = $PT->{'ToLevel'};
  my $PL2 = $PT->{'FromLevel'};
  my $d1 = abs($PL1->{'E'} + $PL1->{'Shift'} - $PL1->{'Er'});
  my $d2 = abs($PL2->{'E'} + $PL1->{'Shift'} - $PL2->{'Er'});
  my ($PLmax,$PLmin,$d_max,$d_min) = ($PL2,$PL1,$d2,$d1);
  if ( $d1 > $d2 ) {
    ($PLmax,$PLmin,$d_max,$d_min) = ($PL1,$PL2,$d1,$d2);
  }
  if ( !$PLmax->{'Precision'} && ($d_max/$Etr1 >= 0.2*$lsd/$lambda1) ) {
    $PLmax->{'Precision'} = 1;
    $NeedTuning = 1;
    $PLmax->{'Shift'} = 0;
    RoundLevel($PLmax);
    my ($lambda1,$u,$lsd) = trans_CalcValueOfObsQuantityRounded($PT);
    $lsd *= 0.5;
    if (!$PLmin->{'Precision'} && (abs($lambda1-$lambda) >= $lsd) && ($d_min/$Etr1 >= 0.2*$lsd/$lambda1)) {
      $PLmin->{'Precision'} = 1;
      $NeedTuning = 1;
      $PLmin->{'Shift'} = 0;
      RoundLevel($PLmin);
    }
  }
} ##trans_tune_SL_levels($)

############################################################################
sub TuneSingleLineLevels {
############################################################################
  # Tune the rounding precision of single-line levels, so that their
  # transitions' Ritz values of the observed quantity exactly match the
  # observed quantity values.

  while ( $NeedTuning ) {
    $NeedTuning = 0;
    for (my $i = 0; $i <= $#SLLines; $i++) {
      my $PT = $SLLines[$i];
      trans_tune_SL_levels($PT);
    }
    if ( $NeedTuning ) {
        # The level shifts have to be re-calculated if tuning is needed,
        # because the tuning changes the level rounding precision,
        # and the shifts are determined by the rounding errors.
      $NeedAdjustment = 1;
      AdjustSingleLineLevels();
    }
  }
}

sub LeftPad($$$) {
  # Left-pad the given string with the given character to the given length.
  my ($S,$PadCh,$Len) = @_;
  my $pad_len = $Len-length($S);
  for (my $i = 1; $i <= $pad_len; $i++) {
    $S = $PadCh . $S;
  }
  return $S;
}

sub RightPad($$$) {
  # Right-pad the given string with the given character to the given length.
  my ($S,$PadCh,$Len) = @_;
  my $pad_len = $Len-length($S);
  for (my $i = 1; $i <= $pad_len; $i++) {
    $S .= $PadCh;
  }
  return $S;
}

sub NotAvail($$) {
  # Return a string with an underscore character _ aligned to a place
  # where the decimal point is printed for decimal numbers. The _ is
  # left-padded and rigt-padded with spaces to the given width.
  my ($width,$places) = @_;
  return LeftPad('_',' ',$width - $places - ($places ? 0 : 1)) .
          RightPad(' ',' ',$places);
}

sub RoundLines() {
  # Properly round the ouput values for each transition.
  my ($Sflags,$PT,$WnCalc,$WnCalcRounded,$UncWnCalc,$WLairCalc,$WLvacCalc,
    $UncWLcalc,$UncWLvacCalc,$UncWLairCalc,$UncWLCalc,$Wobs,$Wunc,$WobsR,$WuncR);

  for (my $i = 1; $i <= $#Lines + 1; $i++) {
    $PT = $Lines[$i-1];
    $WLfactor = $WLfactors{$PT->{'Wunit'}};

    $Sflags = '';
    if ( $PT->{'Flags'} & $FVirtual ) {
      if (!$WriteVirt) {
        next;
      } else {
        $Sflags = 'virt';
      }
    }
    if ( $PT->{'Flags'} & $FExcluded ) {
      $Sflags = 'excl';
    } else {
      $Sflags .= ($Sflags ? ',' : '') . 'P' if ( $PT->{'Flags'} & $FPredict );
      $Sflags .= ($Sflags ? ',' : '') . 'M' if ( $PT->{'Flags'} & $FMask );
      $Sflags .= ($Sflags ? ',' : '') . 'Q' if ( $PT->{'Flags'} & $FQuest );
      $Sflags .= ($Sflags ? ',' : '') . 'b' if ( $PT->{'Flags'} & $FBlend );
      if ( (abs(abs($PT->{'Eobs'}) - $MinAirWn) <= 0.5e-3*$MinAirWn)
           || (abs(abs($PT->{'Eobs'}) - $MaxAirWn) <= 0.5e-3*$MaxAirWn)
      ) {
        if ($PT->{'Flags'} & $Fair) {
           # if observed wavelength is given in air,
           # and it is too close to any end of air range,
           # indicate the air origin of the printed observed wavelength,
           # just to avoid any confusion
          $Sflags .= ($Sflags ? ',' : '') . 'a';
        } else {
          if (!AirWavenumber($PT->{'Eobs'})) {
           # if outside air range, but too close to any end of air range,
           # indicate the vacuum origin of the printed observed wavelength,
           # just to avoid any confusion
            $Sflags .= ($Sflags ? ',' : '') . 'v';
          }
        }
      }
    }
    $PT->{'_flags'} = $Sflags;

    $Wobs = trans_ObsQuantityValue($PT);
    $Wunc = trans_ObsQuantityUnc($PT);
    ($WobsR,$WuncR) = trans_ObsQuantityRounded($PT);

    ($PT->{'_Lobs'},$PT->{'_LuncObs'}) = $RoundInput ? ($WobsR,$WuncR) : ($Wobs,$Wunc);
    ($PT->{'_Eobs'},$PT->{'_WnUncObs'}) = trans_WnObs_rounded($PT);

    ($WnCalc,$UncWnCalc) = trans_Ecalc_from_rounded_levels($PT);
    ($WnCalcRounded,$UncWnCalc) = trans_Ecalc_rounded($PT);
    if ( (abs($WnCalc) < $too_small) && !($PT->{'Flags'} & $FVirtual) ) {
      ($WLairCalc,$WLvacCalc,$UncWLcalc) = ('_','_', '_');
    } else {
      unless ( $SkipLambdas || ($PT->{'Flags'} & $FVirtual) ) {
        ($WLvacCalc,$UncWLvacCalc) = trans_WLvac_calc_rounded($PT);
        $UncWLCalc = $UncWLvacCalc;
        if ( trans_air($PT) ) {
          ($WLairCalc,$UncWLairCalc) = trans_WLair_calc_rounded($PT);
          if ( length($UncWLairCalc) != length($UncWLvacCalc) ) {
              # Force both air and vac wavelength to be rounded to the
              # same number of significant figures, retaining the max.
              # value of their uncertainties
            my $WnUnc = $PT->{'WnUncCalc'};
            my $lambda = 1e8/$WnCalc;
            my $WLunc = ( $Uncert ? abs(($WnUnc/$WnCalc)*$lambda) : 0.001 );
            ($WLvacCalc,$UncWLvacCalc) = round($lambda,$WLunc,$TrRoundThreshold,0);
            $UncWLCalc = $UncWLvacCalc;
            $lambda = Lair(1e8/abs($WnCalc));
            $lambda = -$lambda if ($WnCalc < 0);
            ($WLairCalc,$UncWLairCalc) = round($lambda,$WLunc,$TrRoundThreshold,0);
          }
          $PT->{'_LairCalc'} = $WLairCalc;
          $UncWLCalc = $UncWLairCalc if $UncWLairCalc > $UncWLvacCalc;
        }
        $PT->{'_LvacCalc'} = $WLvacCalc;
        $PT->{'_LuncCalc'} = $UncWLCalc;
      }
      ($PT->{'_WnCalc'},$PT->{'_WnUncCalc'}) = ($WnCalcRounded,$UncWnCalc);

      if ( $Uncert && $PrintUncCorr && !($PT->{'Flags'} & $FVirtual) ) {
          # Format the values of contributions to calculated wavenumber
          # uncertainty due to possible systematic shifts of correlated
          # groups of lines
        $PT->{'_Corr'} = [];
        if (trans_SingleLineLevelFlag($PT) ne '*') {
          for (my $j=0; $j < $MaxCGNo; $j++) {
            ($PT->{'_Corr'}->[$j],$UncWnCalc) =
              round($PT->{'C'}->[$j],$UncWnCalc,$TrRoundThreshold,0);
          }
        } else {
          for (my $j=0; $j < $MaxCGNo; $j++) {
            $PT->{'_Corr'}->[$j] = '';
          }
        }
      }

      unless ( $PT->{'Flags'} & ($FMP | $FVirtual | $FExcluded) ) {
        unless ( $SkipLambdas ) {
          my ($Wc,$WcUnc) = trans_CalcValueOfObsQuantityFromRoundedLevels($PT);
          my $dWunc = $WcUnc;
          $dWunc = $Wunc if ($Wunc > $dWunc);
          ($PT->{'_dW'},$dWunc) = round($Wobs-$Wc,$dWunc,$TrRoundThreshold,0);
        }
        my $dEunc = $PT->{'_WnUncObs'};
        $dEunc = $UncWnCalc if ($UncWnCalc > $dEunc);
        ($PT->{'_dE'},$dEunc) = round($PT->{'Eobs'} - $WnCalc,$dEunc,$TrRoundThreshold,0);
      }
    }
  }
}

sub WriteLines() {
  # Write the output transitions file.
  my $width_E         = $Widths{'widthE'};
  my $width_dW        = $Widths{'width_dW'};
  my $width_dE        = $Widths{'width_dE'};
  my $width_Lobs      = $Widths{'width_Lobs'};
  my $width_Eobs      = $Widths{'width_Eobs'};
  my $width_LuncObs   = $Widths{'width_LuncObs'};
  my $width_WnUncObs  = $Widths{'width_WnUncObs'};
  my $width_LairCalc  = $Widths{'width_LairCalc'};
  my $width_LvacCalc  = $Widths{'width_LvacCalc'};
  my $width_LuncCalc  = $Widths{'width_LuncCalc'};
  my $width_WnCalc    = $Widths{'width_WnCalc'};
  my $width_WnUncCalc = $Widths{'width_WnUncCalc'};
  my $width_flags     = $Widths{'width_flags'};
  my $width_Iobs      = $Widths{'width_Iobs'};
  $width_Iobs         = 5 if $width_Iobs < 5;
  $width_Lobs         = 6 if $width_Lobs < 6;
  $width_LuncObs      = 7 if $width_LuncObs < 7;
  $width_Eobs         = 4 if $width_Eobs < 4;
  $width_WnUncObs     = 7 if $width_WnUncObs < 7;
  $width_LairCalc     = 8 if $width_LairCalc < 8;
  $width_LvacCalc     = 8 if $width_LvacCalc < 8;
  $width_LuncCalc     = 7 if $width_LuncCalc < 7;
  $width_WnCalc       = 5 if $width_WnCalc < 5;
  $width_WnUncCalc    = 7 if $width_WnUncCalc < 7;
  $width_E            = 3 if $width_E < 3;
  $width_dW           = 6 if $width_dW < 6;
  $width_dE           = 6 if $width_dE < 6;
  $width_flags        = 2 if $width_flags < 2;
  $width_flags++;

  my $pl_E         = $Widths{'plE'};
  my $pl_dW        = $Widths{'pl_dW'};
  my $pl_dE        = $Widths{'pl_dE'};
  my $pl_Lobs      = $Widths{'pl_Lobs'};
  my $pl_Eobs      = $Widths{'pl_Eobs'};
  my $pl_LuncObs   = $Widths{'pl_LuncObs'};
  my $pl_WnUncObs  = $Widths{'pl_WnUncObs'};
  my $pl_LairCalc  = $Widths{'pl_LairCalc'};
  my $pl_LvacCalc  = $Widths{'pl_LvacCalc'};
  my $pl_LuncCalc  = $Widths{'pl_LuncCalc'};
  my $pl_WnCalc    = $Widths{'pl_WnCalc'};
  my $pl_WnUncCalc = $Widths{'pl_WnUncCalc'};
  my $Scorr = '';

  open(Tfil1,">$TransOutFilename") or
    Error("Unable to create transitions output file $TransOutFilename");

  # print file header
  my $header = '';
  if ($TabDelimitedOutput) {
    if ($PrintUncCorr) {
      for (my $j=1; $j <= $MaxCGNo; $j++) { $Scorr .= "\tDc$j"; }
    }
    $header .= "Iobs\tW_obs\tuncW_o\twn_o\tuncWnO\t";
    $header .= "W_c_air\tW_c_vac\t" unless $SkipLambdas;
    $header .= "S\t";
    $header .= "uncW_c\t" unless $SkipLambdas;
    $header .= "Wn_c\tuncWnC\t";
    $header .= "dWO-C\t" unless $SkipLambdas;
    $header .= "dEO-C\tL1\tL2\tE1\tE2\tF\tWeight\tUnit$Scorr\n";
  } else {
    if ($PrintUncCorr) {
      for (my $j=1; $j <= $MaxCGNo; $j++) {
        $Scorr .= PrintRound("Dc$j",$width_WnUncCalc,$pl_WnUncCalc);
      }
    }
    $header .=  RightPad('Iobs',   ' ',$width_Iobs + 1) .
                RightPad('W_obs',  ' ',$width_Lobs) .
                RightPad('uncW_o', ' ',$width_LuncObs) .
                RightPad('wn_o',   ' ',$width_Eobs) .
                RightPad('uncWnO', ' ',$width_WnUncObs);
    $header .=  RightPad('W_c_air',' ',$width_LairCalc) .
                RightPad('W_c_vac',' ',$width_LvacCalc) unless $SkipLambdas;
    $header .=  'S ';
    $header .=  RightPad('uncW_c', ' ',$width_LuncCalc) unless $SkipLambdas;
    $header .=  RightPad('Wn_c',   ' ',$width_WnCalc) .
                RightPad('uncWnC', ' ',$width_WnUncCalc);
    $header .=  RightPad('dWO-C',  ' ',$width_dW) unless $SkipLambdas;
    $header .=  RightPad('dEO-C',  ' ',$width_dE) .
                RightPad('L1',     ' ',$LNameWidth) . ' - ' .
                RightPad('L2',     ' ',$LNameWidth + 1) .
                RightPad('E1',     ' ',$width_E) . '- ' .
                RightPad('E2',     ' ',$width_E) .
                RightPad('F',      ' ',$width_flags) .
               'Weight ' .
               'Unit' .
               $Scorr . "\n\n";
  }
  print Tfil1 $header or Error('Error writing to output transitions file.');

  for (my $i = 1; $i <= $#Lines + 1; $i++) {
    my $PT = $Lines[$i-1];
    # The following line added by A. Kramida Aug 25 2017
    next if ( ($PT->{'Flags'} & $FVirtual) && !$WriteVirt);

    my $dW        = $PT->{'_dW'};
    my $dE        = $PT->{'_dE'};
    my $Lobs      = $PT->{'_Lobs'};
    my $Eobs      = $PT->{'_Eobs'};
    my $LuncObs   = $PT->{'_LuncObs'};
    my $WnUncObs  = $PT->{'_WnUncObs'};
    my $LairCalc  = $PT->{'_LairCalc'};
    my $LvacCalc  = $PT->{'_LvacCalc'};
    my $LuncCalc  = $PT->{'_LuncCalc'};
    my $WnCalc    = $PT->{'_WnCalc'};
    my $WnUncCalc = $PT->{'_WnUncCalc'};
    my $flags     = $PT->{'_flags'};
    my $Weight    = sprintf("%4.3f",sqr($PT->{'Weight'}));
    my $unit      = $PT->{'Wunit'};
    $Scorr = '';

    if ($PrintUncCorr) {
      for (my $j=0; $j < $MaxCGNo; $j++) {
        my $Sc = (trans_SingleLineLevelFlag($PT) ne '*') ? $PT->{'_Corr'}->[$j] : '';
        if ( $TabDelimitedOutput ) {
          $Scorr .= "\t$Sc";
        } else {
          $Scorr .= PrintRound($Sc,$width_WnUncCalc,$pl_WnUncCalc);
        }
      }
    }

    my $PL1 = $PT->{'ToLevel'};
    my $PL2 = $PT->{'FromLevel'};

    my $s = '';
    if ($TabDelimitedOutput) {
      my $Slev1 = Trim($PL1->{'Name'});
      my $Slev2 = Trim($PL2->{'Name'});

      $s .= $PT->{'Iobs'} . "\t$Lobs\t$LuncObs\t$Eobs\t$WnUncObs\t";
      $s .= "$LairCalc\t$LvacCalc\t" unless $SkipLambdas;
      $s .= trans_SingleLineLevelFlag($PT) . "\t";
      $s .= "$LuncCalc\t" unless $SkipLambdas;
      $s .= "$WnCalc\t$WnUncCalc\t";
      $s .= "$dW\t" unless $SkipLambdas;
      $s .= "$dE\t$Slev1\t$Slev2\t" . $PL1->{'Er'} . "\t" . $PL2->{'Er'} .
         "\t$flags\t$Weight\t$unit$Scorr\n";
    } else {
      $Lobs      = PrintRound($Lobs,$width_Lobs,$pl_Lobs);
      $Eobs      = PrintRound($Eobs,$width_Eobs,$pl_Eobs);
      $LuncObs   = PrintRound($LuncObs,$width_LuncObs,$pl_LuncObs);
      $WnUncObs  = PrintRound($WnUncObs,$width_WnUncObs,$pl_WnUncObs);
      $LairCalc  = PrintRound($LairCalc,$width_LairCalc,$pl_LairCalc);
      $LvacCalc  = PrintRound($LvacCalc,$width_LvacCalc,$pl_LvacCalc);
      $LuncCalc  = PrintRound($LuncCalc,$width_LuncCalc,$pl_LuncCalc);
      $WnCalc    = PrintRound($WnCalc,$width_WnCalc,$pl_WnCalc);
      $WnUncCalc = PrintRound($WnUncCalc,$width_WnUncCalc,$pl_WnUncCalc);
      $dW        = PrintRound($dW,$width_dW,$pl_dW);
      $dE        = PrintRound($dE,$width_dE,$pl_dE);
      $Scorr     = " $Scorr" if $Scorr;
      $unit      = RightPad($unit,' ',4);

      $s .= sprintf("%-${width_Iobs}s", $PT->{'Iobs'}) .
            $Lobs . $LuncObs . $Eobs . $WnUncObs;
      $s .= $LairCalc . $LvacCalc unless $SkipLambdas;
      $s .= ' ' . trans_SingleLineLevelFlag($PT);
      $s .= $LuncCalc unless $SkipLambdas;
      $s .= $WnCalc .$WnUncCalc;
      $s .= $dW unless $SkipLambdas;
      $s .= $dE .' ' . sprintf("%-${LNameWidth}s", $PL1->{'Name'}) .
            sprintf(" - %-${LNameWidth}s", $PL2->{'Name'}) .
            PrintRound($PL1->{'Er'}, $width_E, $pl_E) . ' -' .
            PrintRound($PL2->{'Er'}, $width_E, $pl_E) .
            sprintf(" %-${width_flags}s", $flags) .
            RightPad($Weight,' ',7) . "$unit$Scorr\n";
    }
    print Tfil1 $s or Error('Error writing to output transitions file.');
  }
  close Tfil1 or Error('Error closing transitions output file.');
}

sub RoundLevel($) {
  # Properly round the output values for the given level.
  my $PL = shift;
  my ($x, $min_unc);

  my $factor = $PL->{'Precision'};
  $factor = "1e$factor";
  ($PL->{'D1r'},$x) = round($PL->{'D1'},$PL->{'D1'},$LevRoundThreshold*$factor,0);
  ($PL->{'Dr'},$x)  = round($PL->{'D'}, $PL->{'D'}, $LevRoundThreshold*$factor,0);
  ($PL->{'D2r'},$x) = round($PL->{'D2'},$PL->{'D2'},$LevRoundThreshold*$factor,0);

   # Choose the minimum of $D, $D1 and $D2 as a criterium for rounding of E
  $min_unc = 'Dr';       # if $PL->{'D'}==0 then choose $PL->{'D1'}
  if (($PL->{$min_unc} == 0) || ($PL->{'D1r'} > 0) && ($PL->{'D1r'} < $PL->{$min_unc})) {
    $min_unc = 'D1r';
  }
  if (($PL->{'D2r'} > 0) && ($PL->{'D2r'} < $PL->{$min_unc})) {
    $min_unc = 'D2r';
  }
  ($PL->{'Er'},$PL->{$min_unc}) = round($PL->{'E'}+$PL->{'Shift'},$PL->{$min_unc},$LevRoundThreshold*$factor,1);
    # now $PL->{$min_unc} accounts for rounding error in $PL->{'E'}

  my $min_uncert = $PL->{$min_unc};
  foreach my $d (@{$PL->{'UncRelBases'}}) {
    $min_uncert = $d if ($d < $min_uncert) && ($d > 0);
  }
  ($min_uncert,$x) = round($min_uncert,$min_uncert,$LevRoundThreshold*$factor,0);
  ($PL->{'D1r'},$x) = round($PL->{'D1'},$min_uncert,$LevRoundThreshold*$factor,0);
  ($PL->{'Dr'},$x)  = round($PL->{'D'}, $min_uncert, $LevRoundThreshold*$factor,0);
  ($PL->{'D2r'},$x) = round($PL->{'D2'},$min_uncert,$LevRoundThreshold*$factor,0);
  ($PL->{'Er'},$x) = round($PL->{'E'}+$PL->{'Shift'},$min_uncert,$LevRoundThreshold*$factor,1);
  $PL->{'min_unc'} = $min_uncert;

    # Upper bound of error due to limited floating-point precision
  ($PL->{'dr'},$x) = round($PL->{'d'},$min_uncert,$LevRoundThreshold,0);
}

sub RoundLevels() {
  # Properly round the output values for each level.
  my $d_max = 0;

  for (my $i = 0; $i <= $#Levels; $i++) {
    my $PL = $Levels[$i];

    RoundLevel($PL);

    if ( $PL->{'d'} > $d_max ) {
      $d_max = $PL->{'d'};
    }

    for (my $k = 0; $k <= $#UncBases; $k++) {
      my $x;
      ($PL->{'UncRelBases'}->[$k],$x) = round($PL->{'UncRelBases'}->[$k],$PL->{'min_unc'},$LevRoundThreshold,0);
    }
  }
  $d_max = $d_max;
} # RoundLevels

sub WriteLevels() {
  # Write the levels output file.
  open(Lfil,">$LevOutFilename") or
    Error("Unable to create levels output file $LevOutFilename");

  my @Comments = ('','fixed by user','fixed by program');
  my $width_E  = $Widths{'widthE'};
  my $width_D  = $Widths{'widthD'};
  my $width_D1 = $Widths{'widthD1'};
  my $width_D2 = $Widths{'widthD2'};
  my $width_d  = $Widths{'widthd'};
  my $pl_E     = $Widths{'plE'};
  my $pl_D     = $Widths{'plD'};
  my $pl_D1    = $Widths{'plD1'};
  my $pl_D2    = $Widths{'plD2'};
  my $pl_d     = $Widths{'pld'};
  my $Nwidth   = $Widths{'LName'};
  $Nwidth   = 12 if $Nwidth < 12;
  $width_E  = 7  if $width_E < 7;
  $width_D  = 3  if $width_D < 3;
  $width_D1 = 3  if $width_D1 < 3;
  $width_D2 = 3  if $width_D2 < 3;
  $width_d  = 2  if $width_D < 2;

  my @basew = ();
  my @basepl = ();
  for ( my $k = 1; $k <= $#UncBases + 1; $k++ ) {
    my $w = $Widths{'Bases'}->[$k-1];
    $w  = 3  if $w < 3;
    $basew[$k-1] = $w;
    my $pl = $Widths{'BasePlaces'}->[$k-1];
    $basepl[$k-1] = $pl;
  }

  my $n = 0;
  my $S = '';
  my $S4 = '';
  if ($#{$Groups}) {
    $S4 = $#{$Groups} + 1;
    $n = length($S4) + 2;
    $width_D2 += ($n + 1);
    #$n;
  }
  my $Scorr = '';
  my $wc = $width_D;
  my $pc = $pl_D;
  if (! $TabDelimitedOutput) {
    if ($MaxCGNo && $PrintUncCorr) {
      $wc = $width_D1 if $wc < $width_D1;
      $wc = $width_D2 if $wc < $width_D2;
      $wc = 4 if $wc < 4;
      $pc = $pl_D1 if $pc < $pl_D1;
      $pc = $pl_D2 if $pc < $pl_D2;
      for ( my $j = 1; $j <= $MaxCGNo; $j++) {
        $Scorr .= RightPad("Dc$j",' ',$wc);
      }
    }
    $Scorr = " $Scorr" if $Scorr;
    print Lfil
      'Designation', RightPad('',' ',$Nwidth-11),
      ' Energy',     RightPad('',' ',$width_E-7),
      ' D1',         RightPad('',' ',$width_D-3),
      ' D2',         RightPad('',' ',$width_D1-3),
      ' D3',         RightPad('',' ',$width_D2-3);
    print Lfil
      ' d',          RightPad('',' ',$width_d-2) if $statistics_size;

    for ( my $k = 1; $k <= $#UncBases + 1; $k++ ) {
      my $w = $basew[$k-1];
      print Lfil ($k == 1 ? ' ' : '') . "B$k", RightPad('',' ',$w-2);
    }
    print Lfil  $Scorr, "    N_lines         Comments\n\n";
  } else {
    print Lfil  join("\t",(
      'Designation',
      'Energy',
      'D1',
      'D2',
      'D3'));
    print Lfil "\t" .'d' if $statistics_size;

    for ( my $k = 1; $k <= $#UncBases + 1; $k++ ) {
      print Lfil "\tB$k";
    }
    if ($PrintUncCorr) {
      for ( my $j = 1; $j <= $MaxCGNo; $j++) {
        $Scorr .= "\tDc$j";
      }
    }
    print Lfil  $Scorr, "\tN_lines\tComments\n";
  }
  for (my $i = 1; $i <= $#Levels + 1; $i++) {
    my $PL = $Levels[$i-1];
         # Form 'NLines' description
    $S = '    0';
    if ($PL->{'NT1'}) {
      $S = sprintf("%5d",$PL->{'NT2'});
      $S .= ', '.$PL->{'NB'}.'B' if ($PL->{'NB'});
      $S .= ', '.$PL->{'NQ'}.'Q' if ($PL->{'NQ'});
      $S .= ', '.$PL->{'NL'}.'L' if ($PL->{'NL'});
      $S .= ', '.$PL->{'ND'}.'D' if ($PL->{'ND'});
    }
    my $base_unc = '';
    for ( my $k = 1; $k <= $#UncBases + 1; $k++ ) {
      my $unc = $PL->{'UncRelBases'}->[$k-1];
      if ( $unc == '_' ) {
        my $BL = $UncBases[$k-1];
        $unc = "B$k" if ($PL eq $BL);
      }
      if ($TabDelimitedOutput) {
        $S = Trim($S);
        $base_unc .= "\t$unc";
      } else {
        $base_unc .= ($k == 1 ? ' ' : '') . RightPad($unc,' ',$basew[$k-1]);
      }
    }

    $S4 = ($PL->{'GNo'} > 1) ? '['.$PL->{'GNo'}.']' : '' ;
    $S4 = RightPad($S4,' ',$n) unless $TabDelimitedOutput;
    my $dd3 = $PL->{'D'};
    my $E = $TabDelimitedOutput ? $PL->{'Er'} : PrintRound($PL->{'Er'},$width_E,$pl_E);
    $Scorr = '';
    my $Dc = 0.0;

    # Determine presentation parameters $dd4, $p4 and $w4 for corr. unc.
    my $dd4 = $PL->{'D1'};
    my $p4 = num_dec_places($PL->{'D1r'});
    if ($PL->{'D'} < $dd4) {
      $dd4 = $PL->{'D'};
      my $p = num_dec_places($PL->{'Dr'});
      $p4 = $p if ($p > $p4);
    }
    if (($PL->{'D2'} > $too_small) && ($PL->{'D2'} < $dd4)) {
      $dd4 = $PL->{'D2'} ;
      my $p = num_dec_places($PL->{'D2r'});
      $p4 = $p if ($p > $p4);
    }
    if (!$PL->{'Fixed'}) {
      for (my $j = 0; $j < $MaxCGNo; $j++) {
        my $FFc = $PL->{'C'}->[$j];
        if ($PrintUncCorr) {
          my ($Sc,$x) = round($FFc,$dd4,$LevRoundThreshold,0);
          if ( $TabDelimitedOutput ) {
            $Scorr .=  "$Sc\t";
          } else {
            $Scorr .=  PrintRound($Sc,$wc,$p4);
          }
        }
        $Dc += sqr($FFc);
      }

      if ($PL->{'NT2'} > 1) {
        # If not a single-line level, add estimated upper bound
        # of possible systematic error to "absolute uncertainty" estimates
        $PL->{'D1'} = sqrt(sqr($PL->{'D1'}) + $Dc);
        my $x;
        ($PL->{'D1r'},$x) = round($PL->{'D1'},$PL->{'min_unc'},$LevRoundThreshold,0);
        if ($PL->{'GNo'}>1) {
          $PL->{'D2'} = sqrt(sqr($PL->{'D2'}) + $Dc);
          ($PL->{'D2r'},$x) = round($PL->{'D2'},$PL->{'min_unc'},$LevRoundThreshold,0);
        }
      }
    } else {
      if ($PrintUncCorr) {
        if ($TabDelimitedOutput) {
          $Scorr = RightPad('',"\t",$MaxCGNo);
        } else {
          $Scorr = RightPad('',' ',$wc*$MaxCGNo);
        }
      }
    }

    my $S2 = '';
    my $Slev = '';
    if ($PL->{'D2'} > $too_small) {
      $S2 = PrintRound($PL->{'D2r'},$width_D2-$n,$pl_D2) . $S4;
    } else {
      if ($TabDelimitedOutput) {
        $S2  =  '_' . $S4;
      } else {
        $S2 = LeftPad('_',' ', $width_D2 - $pl_D2 + 1 - ($pl_D2 <= 0 ? 1 : 0) - $n) . $S4;
      }
    }

    if ($TabDelimitedOutput) {
      $Slev = Trim($PL->{'Name'});
      $S = Trim($S);
      $S2 .= "\t"  . $PL->{'dr'} if $statistics_size;
      print Lfil "$Slev\t" ,
                 "$E\t",
                 $PL->{'Dr'},"\t",
                 $PL->{'D1r'},"\t",
                 "$S2$base_unc\t$Scorr$S\t",
                 $Comments[$PL->{'Fixed'}],"\n";
    } else {
      my $name_pad = RightPad('',' ',$Nwidth-length($PL->{'Name'}));
      $S2 .= PrintRound($PL->{'dr'},$width_d,$pl_d) if $statistics_size;
      $Scorr = " $Scorr" if $Scorr;
      print Lfil $PL->{'Name'}, $name_pad,
        $E,
        PrintRound($PL->{'Dr'},$width_D,$pl_D),
        PrintRound($PL->{'D1r'},$width_D1,$pl_D1),
        $S2,
        $base_unc,
        $Scorr,
        $S,
        RightPad('',' ',20-length($S)),
        $Comments[$PL->{'Fixed'}],"\n";
    }
  }
  close(Lfil) or Error("Error closing levels output file.");
} # WriteLevels

sub AdjustWeights($$) {
    # Calculate correct number of transitions defining each level
    # and adjust values of dispersion: exclude transitions that are
    # the only one connection for upper or lower level.
  my ($G1,$GNo) = @_;
  my ($i,$ii);
  my $ii = 1;
  do {
    $ii = 0;
    my $lines_count = $#GroupLines + 1;
    for (my $i=1; $i <= $lines_count; $i++) {
      my $PT = $GroupLines[$i-1];
      my $From = $PT->{'FromLevel'};
      my $To = $PT->{'ToLevel'};
      next if ( ($PT->{'Flags'} & $FMP) || ($PT->{'WnUnc'} < $too_small)
        || ($PT->{'Flags'} & $FVirtual) );
      next if ( $To->{'Fixed'} && ($To->{'GNo'}!=$GNo) );
      next if ( $From->{'Fixed'} && ($From->{'GNo'}!=$GNo) );

      if ( !($PT->{'Flags'} & $FExcludedTo)
        && ($From->{'NT1'} + $From->{'NV'} == 1)
        && ($To->{'NT1'} + $To->{'NV'} > 1) )
      {
        $ii = 1;
          # decrease number of transitions
        $To->{'NT1'}--;
        $PT->{'Flags'} |= $FExcludedTo;
        my $dE = $PT->{'WnUnc'};
        my $dE1 = $PT->{'Eobs'} - trans_Ecalc($PT);
          # exclude the weight of this transition from the weights sum
          # for the level $To
        my $FF = 1/sqr($dE);
        my $dd = $FF*(1+sqr($dE1/$dE));
        $To->{'D'} -= $dd;
        $To->{'NT'} -= $FF;

        next if ( $To->{'Fixed'} && ($To->{'GNo'} != $GNo) );

        # decrease number of doubly-assigned lines if the line is doubly-assigned
        $To->{'ND'}-- if trans_DoublyAssigned($PT);
        # decrease number of questionable lines if the line is questionable
        if ($PT->{'Flags'} & $FQuest) {
          $To->{'NQ'}-- ;
        } elsif ($PT->{'Flags'} & $FBlend) {
            # decrease number of blends if the line is blended
          $To->{'NB'}--;
        } else {
            # decrease number of lines with high deviations
          $To->{'NL'}-- if (abs($dE1) > 1.2*abs($dE));
        }
      }
      if ( !($PT->{'Flags'} & $FExcludedFrom)
        && ($To->{'NT1'} + $To->{'NV'} == 1)
        && ($From->{'NT1'} + $From->{'NV'} > 1)
      ) {
        $ii = 1;
          # decrease number of transitions
        $From->{'NT1'}--;
        $PT->{'Flags'} |= $FExcludedFrom;
        my $dE  = $PT->{'WnUnc'};
        my $dE1 = $PT->{'Eobs'} - trans_Ecalc($PT);
          # exclude the weight of this transition from the weights sum
          #  for the level $PT->.FromLevel
        my $FF = 1/sqr($dE);
        my $dd = $FF*(1+sqr($dE1/$dE));
        $From->{'D'} -= $dd;
        $From->{'NT'} -= $FF;
        next if ($From->{'Fixed'} && ($From->{'GNo'} != $GNo) );

        # decrease number of doubly-assigned lines if the line is doubly-assigned
        $From->{'ND'}-- if (trans_DoublyAssigned($PT));

        # decrease number of questionable lines if the line is questionable
        if ($PT->{'Flags'} & $FQuest) {
          $From->{'NQ'}--;
        } elsif ($PT->{'Flags'} & $FBlend) {
            # decrease number of blends if the line is blended
          $From->{'NB'}--;
        } else {
            # decrease number of lines with high deviations
          $From->{'NL'}-- if (abs($dE1)>1.2*abs($dE));
        }
      }
    }
  } while ($ii);
} # AdjustWeights

sub AdjustDispersions ($$) {
     # AdjustDispersions:
     # find final values of Radziemski dispersion ({'D'}, called D1 in the output file).
  my ($G1,$GNo) = @_;
  my $max_d = 0;
  foreach my $PL (@{$G1}) {
    if (!$PL->{'Fixed'} || ($PL->{'GNo'}==$GNo)) {
      if (($PL->{'D'} > $too_small) && ($PL->{'NT'} > $too_small)) {
        $PL->{'D'} = sqrt(($PL->{'D'}/$PL->{'NT'})/$PL->{'NT'});
      } else {
        $PL->{'D'} = 0;
      }
    }

    my $d = ($PL->{'D'} > 0) ? $PL->{'d'}/$PL->{'D'} : 0;
    $max_d = $d if $max_d < $d;

    if (!$PL->{'Fixed'}) {
      $PL->{'NT2'} = $PL->{'NT1'};
    } else {
      if ($PL->{'GNo'}==$GNo) {
        $PL->{'NT2'} = $PL->{'NT1'};
      }
    }
  }
  $max_d = $max_d;
}

sub SortLevels() {
   # Sort levels by energy
  my ($sw,$PL,$PL1,$E1,$E);
  my $sw = 1;
  while ( $sw ) {
    $sw = 0;
    for (my $i = 1; $i <= $#Levels + 1; $i++) {
      for (my $k = $i+1; $k <= $#Levels + 1; $k++) {
        $PL = $Levels[$i-1];
        $PL1 = $Levels[$k-1];
        $E1 = $PL1->{'Er'};
        $E  = $PL->{'Er'};
        if ( ($E1 < $E) or ($E1 == $E) && ($PL1->{'E'} < $PL->{'E'}) ) {
          $Levels[$i-1] = $PL1;
          $Levels[$k-1] = $PL;
          $sw = 1;
        }
      }
    }
  }
}

############################################################################
sub FindRSS($$) {   #10/20/2009 3:03PM
############################################################################
  my ($G, $GNo) = @_;
  my ($deg,$Nlev,$Nlin) = (0,0,0);
  my $RSS = 0;

  foreach my $Lev (@{$G}) {
    next if ($Lev->{'GNo'} != $GNo);
    next if ($Lev eq $Ground);        # The ground level is not counted
    $Nlev++;
  }

  foreach my $PT (@GroupLines) {
    next if ($PT->{'Flags'} & $FMP);
    next if ($PT->{'Flags'} & $FExcluded);
    next if ($PT->{'FromLevel'}->{'GNo'} != $GNo);
    next if ($PT->{'ToLevel'}->{'GNo'} != $GNo);
    $RSS += sqr((trans_Ecalc($PT) - $PT->{'Eobs'})/$PT->{'WnUnc'});
    $Nlin++;
  }
  $deg = ($Nlev ? $Nlin - $Nlev : 0);
  $deg = 0 if $deg < 0;
  $RSS = ($deg ? sprintf("%4.2f",$RSS/$deg) : 0);
  return ($RSS,$deg);
} ##FindRSS($$)

############################################################################
sub FindConditionNumber($$) {  #11/19/2009 2:46PM
############################################################################
  my ($A,$B) = @_;
  my $NA = $#{$A->[0]} + 1;
  my $norm_Fro_A = 0;
  my $norm_Fro_B = 0;
  my $norm_row_max_A = 0;
  my $norm_row_max_B = 0;
    # Find Frobenius norms and max-row norms of A and A^(-1) = B
  for ( my $i = 0; $i < $NA; $i++ ) {
    my $row_sum_A = 0;
    my $row_sum_B = 0;
    for ( my $j = 0; $j < $NA; $j++ ) {
      $row_sum_A += abs($A->[$i]->[$j]);
      $row_sum_B += abs($B->[$i]->[$j]);
      $norm_Fro_A += sqr($A->[$i]->[$j]);
      $norm_Fro_B += sqr($B->[$i]->[$j]);
    }
    $norm_row_max_A = $row_sum_A if $row_sum_A > $norm_row_max_A;
    $norm_row_max_B = $row_sum_B if $row_sum_B > $norm_row_max_B;
  }
  $norm_Fro_A = sqrt($norm_Fro_A);
  $norm_Fro_B = sqrt($norm_Fro_B);

    # Find max-col norms of A and A^(-1) = B
  my $norm_col_max_A = 0;
  my $norm_col_max_B = 0;
  for ( my $j = 0; $j < $NA; $j++ ) {
    my $col_sum_A = 0;
    my $col_sum_B = 0;
    for ( my $i = 0; $i < $NA; $i++ ) {
      $col_sum_A += abs($A->[$i]->[$j]);
      $col_sum_B += abs($B->[$i]->[$j]);
    }
    $norm_col_max_A = $col_sum_A if $col_sum_A > $norm_col_max_A;
    $norm_col_max_B = $col_sum_B if $col_sum_B > $norm_col_max_B;
  }

  return ($norm_row_max_A*$norm_row_max_B);
} ##FindConditionNumber($$)

############################################################################
sub FindRelCondNumber($$) {  #12/06/2009 9:30AM
############################################################################
  my ($A,$B) = @_;    # Here B must be equal to A^-1
  my $NA = $#{$A->[0]} + 1;
  my $row_sum = 0;
  my $norm_row_max = 0;
  for ( my $i = 0; $i < $NA; $i++ ) {
    $row_sum = 0;
    my $BAij = 0;
    for ( my $j = 0; $j < $NA; $j++ ) {
      for ( my $k = 0; $k < $NA; $k++ ) {
        $BAij += abs($B->[$i]->[$k])*abs($A->[$k]->[$j]);
      }
      $row_sum += $BAij;
    }
    $norm_row_max = $row_sum if $norm_row_max < $row_sum;
  }
  return $norm_row_max;
} ##FindRelCondNumber($$)

############################################################################
#                       Main program
############################################################################
$Groups = [];
@Lines = ();
@VirtLines = ();
@FixLevs = ();
@Levels = ();

my $time_in = time();

print "Reading parameters...";
&ReadParams(shift);
print "Done.\n";
print $use_Cholesky ? "\nUsing Cholesky decomposition to solve the matrix equation.\n\n" :
  "\nUsing LU decomposition to solve the matrix equation.\n\n";

   # Read transitions wavelengths and designations from text file
print "Reading lines...";
&ReadLines();
my $t0 = time();
print "Done.\n";

print "Filling groups...";
&FillCGroups();
print "Done.\n";

print "Counting assigned lines...";
&FindNumAssign();
print "Done.\n";

print "Reading fixed levels...";
&ReadFixLevs();
print "Done.\n";
my $t1 = time();

for (my $GNo = 1; $GNo <= $#{$Groups} + 1; $GNo++) {
  @GroupLines = ();
  my $G1 = $Groups->[$GNo-1];
  print "\nGroup $GNo of ",$#{$Groups}+1,': ',$#{$G1}+1," levels.\n";
  print "Pass 1: finding energies...\n";

          # Find energies
  FormMatrix($G1,$GNo,1);
  FindLevels($G1,$GNo,1);
  $t1 = time();
  print "Done. (spent ", $t1 - $t0, " sec)\n";
  print "Number of substituted variables = $num_subst\n";
  print "Pass 2: finding dispersions...\n";

        # Find dispersions
  FormMatrix($G1,$GNo,2);
  FindLevels($G1,$GNo,2);

          # Calculate number of lines with high deviations
          #  (but not blended or questionable)
          #  for each level
  for (my $j = 0; $j <= $#GroupLines; $j++) {
    my $PT = $GroupLines[$j];
    next if ($PT->{'Flags'} & $FMP);

  # increase number of doubly-assigned lines if the line is doubly-assigned
    if (trans_DoublyAssigned($PT)) {
      $PT->{'ToLevel'}->{'ND'}++;
      $PT->{'FromLevel'}->{'ND'}++;
    }
  # increase number of large-deviation lines if the line deviates too much
    if ( !($PT->{'Flags'} & ($FBlend | $FQuest)) &&
      (abs($PT->{'Eobs'} - trans_Ecalc($PT)) > 1.2*abs($PT->{'WnUnc'}))
    ) {
      $PT->{'FromLevel'}->{'NL'}++;
      $PT->{'ToLevel'}->{'NL'}++;
    }
  }
  $t0 = time();
  print "Done. (spent ", $t0 - $t1, " sec)\n";

  print "Adjusting weights...";
  AdjustWeights($G1,$GNo);
  print "Done.\n";

  print "Adjusting dispersions...";
  AdjustDispersions($G1,$GNo);
  print "Done.\n";

  my ($RSS,$degrees_of_freedom) = FindRSS($G1,$GNo);
  print "\nRSS/degrees_of_freedom = $RSS ($degrees_of_freedom degrees_of_freedom)\n\n";

  if ($Uncert) {
    print "Finding Ritz uncertainties for lines...\n";
    FindUncert($G1,$GNo);
    my $t = time();
    print "Done. (spent ", $t - $t0, " sec)\n";
    $t0 = $t;
  }

  if ( $max_num_error_frac > 0.5 ) {
    FormMatrix($G1,$GNo,1);
    my ($rms,$big_dev) = (0,0);
    my $A_scaled = scale_rows($A);
    ($B,$rms,$big_dev) = InvertMatrix($A_scaled,$InvChk,0);
    $cond = &FindConditionNumber($A_scaled,$B);
    #$rel_cond = &FindRelCondNumber($A_scaled,$B);
    $normwise_error_bound = $eps*$cond;
    $normwise_error_bound /= (1-$eps*$cond/2) if ($eps*$cond/2 < 1);
    print sprintf("Max. ratio of estimated numerical errors to level uncertainties = %.1e\n",$max_num_error_frac);
    print "Warning: numerical errors due to finite floating-point precision may be significant in some levels.\n";
    print "Check column 'd' in the levels output file. The next four lines contain relevant diagnostics.\n";
    print sprintf("Matrix condition number         = %.1e\n",$cond);
    print sprintf("Machine precision               = %.1e\n",$eps);
    print sprintf("Normwise relative error bound   = %.1e\n",$normwise_error_bound);
    #print sprintf("Relative condition number      = %6.1e\n",$rel_cond);
    #print sprintf("Iteration improvement criterion = %6.1e\n\n",$iter_crit);
  }
}

print "\nRounding levels...";
RoundLevels();
print "Done.\n";

print "Adjusting single-line levels...";
AdjustSingleLineLevels();
print "Done.\n";

if ( $TuneLevels ) {
  print "Tuning single-line levels...";
  TuneSingleLineLevels();
  print "Done.\n";
}

print "Rounding lines...";
RoundLines();
print "Done.\n";

unless ( $TabDelimitedOutput ) {
  print "Finding column widths...";
  FindWidths;
  print "Done.\n";
}

print "Sorting levels...";
SortLevels;
print "Done.\n";

print "Writing levels...";
WriteLevels;
print "Done.\n";

print "Writing lines...";
WriteLines;
print "Done.\n";

print("Ok.\n");
$t1 = time() - $time_in;
print "Total time spent: $t1 sec\n";

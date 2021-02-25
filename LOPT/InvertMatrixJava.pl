use strict;

############################################################################
sub transpose($) {   #01/18/2011 3:36PM
############################################################################
  my $A = shift;
  my $rows = $#{$A} + 1;
  my $cols = $#{$A->[0]} + 1;
  my $B = [];
  for ( my $j = 0; $j < $cols; $j++ ) {
    $B->[$j] = [];
  }
  for ( my $i = 0; $i < $rows; $i++ ) {
    for ( my $j = 0; $j < $cols; $j++ ) {
      $B->[$j]->[$i] = $A->[$i]->[$j];
    }
  }
  return $B;
} ##transpose($)

############################################################################
sub LUDecomposition($) {   #12/08/2009 6:53PM
############################################################################
  # Calculate an LU decomposition of a matrix.
  #@param $A   square matrix n x n
  #@return     LU(r(:),:)
  #
  # Adapted from JAMA package, J. Hicklin, C. Moler, P. Webb, R.F. Boisvert,
  # B. Miller, R. Pozo, K. Remington, JAMA: A Java Matrix Package,
  # online at http://math.nist.gov/javanumerics/jama/ (2005),

  my $A = shift;
  my $N = $#{$A} + 1;

  return if ($N<=0);

  open(M_OUT, ">matrix_in.txt") or die ("Could not create file matrix_in.txt");

  print M_OUT &MatrixToString($A,' ');
  close M_OUT or die "Error closing file matrix_in.txt";

  my $java_output = `java -jar InvertMatrix.jar LU 0 0`;
  print $java_output;
  return &LUDecomposition1($A) if ( $java_output =~ /error/i );

  open(M_IN, "<matrix_out.txt") or die ("Could not open file matrix_out.txt");

  my $result = [];
  my $s;
  my $i = -2;
  while ( ($s = <M_IN>) && (defined $s) ) {
    $i++;
    last if ($i == $N);
    next if $i < 0; # skip the first line
    chomp $s;
    next unless $s;
    $s =~ s/^\s+|\s+$//;
    my @row = split(/ +/, $s);
    $result->[$i] = \@row;
  }
  if ( $i != $N ) {
    die "Matrix size after LU decomposion $i does not match the input matrix size $N\n";
  }
  my $pivot = [];
  while ( ($s = <M_IN>) && (defined $s) ) {
    chomp $s;
    $s =~ s/^\s+|\s+$//;
    next unless $s;
    my @row = split(/ +/, $s);
    $i = $#row + 1;
    if ( $i != $N ) {
      die "Pivot size after LU decomposion $i does not match the input matrix size $N\n";
    }
    for ( my $j = 0; $j < $N; $j++ ) {
      my $p = int($row[$j] + 0.1);
      $row[$j] = $p;
    }
    $pivot = \@row;
    last;
  }
  close(M_IN);
  return ($result,$pivot);
} ##LUDecomposition($)

sub LUDecomposition1 ($) {
  # Calculate an LU decomposition of a matrix.
  #@param $A   square matrix n x n
  #@return     LU(r(:),:)

  # Adapted from disabled experimental code in JAMA package,
  # J. Hicklin, C. Moler, P. Webb, R.F. Boisvert,
  # B. Miller, R. Pozo, K. Remington, JAMA: A Java Matrix Package,
  # online at http://math.nist.gov/javanumerics/jama/ (2005),
  # and optimized for Perl.

  my $A = shift;
      # Initialize.
  my @LU = ();
  my $n = $#{$A} + 1; # dimension of square matrix $A
  my @piv = ();
  for (my $i = 0; $i < $n; $i++) {
    $piv[$i] = $i;
    $LU[$i] = [];
    push(@{$LU[$i]},@{$A->[$i]});
  }
  #my $pivsign = 1;

  # Main loop.
  for (my $k = 0; $k < $n; $k++) {
    # Find pivot.
    my $p = $k;
    for (my $i = $k+1; $i < $n; $i++) {
      $p = $i if (abs($LU[$i]->[$k]) > abs($LU[$p]->[$k]));
    }
    # Exchange if necessary.
    if ($p != $k) {
      my $t = $LU[$p]; $LU[$p] = $LU[$k]; $LU[$k] = $t;
      $t = $piv[$p]; $piv[$p] = $piv[$k]; $piv[$k] = $t;
      #$pivsign = -$pivsign;
    }
    # Compute multipliers and eliminate k-th column.
    my $LUk = $LU[$k];
    if ($LUk->[$k]) {
      my $LUkk = 1/$LUk->[$k];
      for (my $i = $k+1; $i < $n; $i++) {
        my $LUi = $LU[$i];
        # Take advantage of sparsity of $A: skip the column cycle
        # if $LUi->[$k] is zero.
        if ( $LUi->[$k] ) {
          my $LUik = $LUi->[$k] *= $LUkk;
          for (my $j = $k+1; $j < $n; $j++) {
            $LUi->[$j] -= $LUk->[$j]*$LUik;
          }
        }
      }
    } else {
      # Error: matrix is not invertible. Return the index of error row
      return ($k);
    }
  }
  return (\@LU,\@piv);
} ##LUDecomposition1($)

############################################################################
sub LUSolve($$$) {  #12/09/2009 2:51PM
############################################################################
 #     Solve A*X = B  where A = LU
 #  param   $LU  LU decomposition of a square matrix n x n
 #  param   $B   A Matrix with as many rows as A and any number of columns.
 #  @return     X so that L*U*X = B(piv,:)

  # Adapted from JAMA package, J. Hicklin, C. Moler, P. Webb, R.F. Boisvert,
  # B. Miller, R. Pozo, K. Remington, JAMA: A Java Matrix Package,
  # online at http://math.nist.gov/javanumerics/jama/ (2005),
  # and optimized for Perl.

  my ($LU,$B,$piv) = @_;
  my $n = $#{$LU} + 1; # number of rows of square matrix $LU and of $B
  my $nx = $#{$B->[0]} + 1; # number of columns of $B
  return if (($n<=0) || ($nx<=0));

  open(M_OUT, ">matrix_in.txt") or die ("Could not create file matrix_in.txt");

  print M_OUT &MatrixToString($LU,' ');
  print M_OUT &MatrixToString(&transpose($B),' ');
  print M_OUT join(' ',@{$piv}) . "\n";
  close M_OUT or die "Error closing file matrix_in.txt";

  my $java_output = `java -jar InvertMatrix.jar LUS 0 0`;
  print $java_output;
  return '' if ( $java_output =~ /error/i );

  open(M_IN, "<matrix_out.txt") or die ("Could not open file matrix_out.txt");

  my $result = [];
  my $s;
  my $i = -2;
  while ( ($s = <M_IN>) && (defined $s) ) {
    $i++;
    next if $i < 0; # skip the first line
    chomp $s;
    next unless $s;
    $s =~ s/^\s+|\s+$//;
    my @row = split(/ +/, $s);
    $result->[$i] = \@row;
  }
  if ( $i != $n ) {
    die "Matrix size after LUSolve $i does not match the input matrix size $n\n";
  }
  close(M_IN);
  return $result;
} ##LUSolve($$)

############################################################################
sub CholeskyDecomposition($) {
#    Cholesky algorithm for symmetric and positive definite matrix.
#   @param  A   Square, symmetric matrix.
#   @return     Structure to access L and isspd flag.
############################################################################
  my $A = shift;
      # Initialize.
  my $N = $#{$A} + 1; # dimension of square matrix $A
  return if ($N<=0);

  open(M_OUT, ">matrix_in.txt") or die ("Could not create file matrix_in.txt");

  print M_OUT &MatrixToString($A,' ');
  close M_OUT or die "Error closing file matrix_in.txt";

  my $java_output = `java -jar InvertMatrix.jar CH 0 0`;
  print $java_output;
  return &CholeskyDecomposition1($A) if ( $java_output =~ /error/i );

  open(M_IN, "<matrix_out.txt") or die ("Could not open file matrix_out.txt");

  my $result = [];
  my $s;
  my $i = -2;
  while ( ($s = <M_IN>) && (defined $s) ) {
    $i++;
    next if $i < 0; # skip the first line
    chomp $s;
    next unless $s;
    $s =~ s/^\s+|\s+$//;
    my @row = split(/ +/, $s);
    $result->[$i] = \@row;
  }
  if ( $i != $N ) {
    die "Matrix size after Cholesky decomposion $i does not match the input matrix size $N\n";
  }
  close(M_IN);
  return $result;
}

############################################################################
sub CholeskyDecomposition1($) {
#    Cholesky algorithm for symmetric and positive definite matrix.
#   @param  A   Square, symmetric matrix.
#   @return     Structure to access L and isspd flag.
############################################################################
  my $A = shift;
      # Initialize.
  my $n = $#{$A} + 1; # dimension of square matrix $A
  my $L = [];
  for ( my $i=0; $i < $n; $i++ ) {
    $L->[$i] = [];
    for ( my $j=0; $j < $n; $j++ ) {
      $L->[$i]->[$j] = 0;
    }
  }
  my $isspd = 1;
  my $asymmetric = 0;
  my $error_row = -1;

      # Main loop.
  for (my $j = 0; $j < $n; $j++) {
    my $Lrowj = $L->[$j];

    my %Lrowj_map = ();
    for (my $k = 0; $k < $j; $k++) {
      if ($Lrowj->[$k]) {
        $Lrowj_map{$k} = 1;
      }
    }

    my $Arowj = $A->[$j];
    my $d = 0.0;
    for (my $k = 0; $k < $j; $k++) {
      my $Lrowk = $L->[$k];
      my $Lkk = $Lrowk->[$k];
      my $s = 0.0;
      foreach my $i (keys %Lrowj_map) {
        $s += $Lrowk->[$i]*$Lrowj->[$i] if $i < $k;
      }
      $Lrowj->[$k] = $s = ($Arowj->[$k] - $s)/$Lkk;
      $Lrowj_map{$k} = 1 if $s;
      $d += $s*$s;
      if ( $A->[$k]->[$j] != $Arowj->[$k] ) {
        $isspd = 0;
        $asymmetric = 1;
        $error_row = -1;
        last;
      }
    }
    $d = $Arowj->[$j] - $d;
    if ( $d <= 0 ) {
      $isspd = 0;
      $error_row = $j unless $asymmetric;
      last;
    }
    $Lrowj->[$j] = sqrt($d);
  }
  return ($isspd) ? $L : $error_row;
}

############################################################################
sub SolveCholesky($$) { #01/28/2010 10:51AM
############################################################################
  my ($L,$B) = @_;
  my $n = $#{$L} + 1; # dimension of square matrix $L -- must be the same as row dimension of $B
  my $nx = $#{$B->[0]} + 1; # column dimension of $B
  return if (($n<=0) || ($nx<=0));

  open(M_OUT, ">matrix_in.txt") or die ("Could not create file matrix_in.txt");

  print M_OUT &MatrixToString($L,' ');
  print M_OUT &MatrixToString(&transpose($B),' ');
  close M_OUT or die "Error closing file matrix_in.txt";

  my $java_output = `java -jar InvertMatrix.jar CHS 0 0`;
  print $java_output;
  return '' if ( $java_output =~ /error/i );

  open(M_IN, "<matrix_out.txt") or die ("Could not open file matrix_out.txt");

  my $result = [];
  my $s;
  my $i = -2;
  while ( ($s = <M_IN>) && (defined $s) ) {
    $i++;
    next if $i < 0; # skip the first line
    chomp $s;
    next unless $s;
    $s =~ s/^\s+|\s+$//;
    my @row = split(/ +/, $s);
    $result->[$i] = \@row;
  }
  if ( $i != $n ) {
    die "Matrix size after SolveCholesky $i does not match the input matrix size $n\n";
  }
  close(M_IN);
  return $result;
} ##SolveCholesky

############################################################################
sub InvertMatrix($$;$) {   #12/10/2009 11:22AM
############################################################################
 #  Invert a matrix A
 #  param    $A       square matrix n x n
 #  param    $InvChk  flag (0 or 1 if invershion check is needed)
 #  @return  A^-1 if A is invertible, $error_row otherwise
 #  @return  $RMS     root-mean-square deviation of A*A^-1 - I from zero if $InvChk==1 and matrix is invertible
 #  @return  $BIG     largest deviation of A*A^-1 - I from zero if $InvChk==1 and matrix is invertible

  my ($A, $InvChk, $UseCholesky) = @_;
  $InvChk += 0;
  $UseCholesky += 0;
  my $N = $#{$A} + 1;
  return if ($N<=0);

  open(M_OUT, ">matrix_in.txt") or die ("Could not create file matrix_in.txt");

  print M_OUT &MatrixToString($A,' ');
  close M_OUT or die "Error closing file matrix_in.txt";

  my $java_output = `java -jar InvertMatrix.jar INV $InvChk $UseCholesky`;
  if ( $java_output =~ /error/i ) {
    print $java_output;
    return '';
  }

  my ($RMS,$BIG) = (0,0);
  if ( $java_output =~ /RMS = ([0-9.E+-]+) +Largest deviation = ([0-9.E+-]+)/mi ) {
    ($RMS,$BIG) = ($1, $2);
  }

  open(M_IN, "<matrix_out.txt") or die ("Could not open file matrix_out.txt");

  my $result = [];
  my $s;
  my $i = -2;
  while ( ($s = <M_IN>) && (defined $s) ) {
    $i++;
    next if $i < 0; # skip the first line
    chomp $s;
    next unless $s;
    $s =~ s/^\s+|\s+$//;
    my @row = split(/ +/, $s);
    $result->[$i] = \@row;
  }
  if ( $i != $N ) {
    die "Matrix size after inversion $i does not match the size before inversion $N\n";
  }
  return ($result,$RMS,$BIG);
} ##InvertMatrix($$;$)

############################################################################
sub MatrixToString($$) {  #10/14/2009 12:08PM
############################################################################
  my ($A, $delim) = @_;
  my $N = $#{$A} + 1;
  my $s = '';
  for (my $row = 0; $row < $N; $row++) {
    $s .= join($delim,@{$A->[$row]}) . "\n";
  }
  return $s;
} ##MatrixToString($)


1;

# Testing code
#open(M_IN, "<matrix_in.txt") or die ("Could not open file matrix_in.txt");
#my $A = [];
#my $s;
#my $i = -1;
#while ( ($s = <M_IN>) && (defined $s) ) {
# $i++;
#  chomp $s;
#  next unless $s;
#  $s =~ s/^\s+|\s+$//;
#  my @row = split(/ +/, $s);
#  $A->[$i] = \@row;
#}
#my $B = InvertMatrix($A,1);

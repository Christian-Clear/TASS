use strict;

############################################################################
sub LUDecomposition($) {   #12/08/2009 6:53PM
############################################################################
  # Calculate an LU decomposition of a matrix.
  #@param $A   square matrix n x n
  #@return     LU(r(:),:)
  #
  # This method is not used in the current version, as it turned out to be
  # twice as slow as LUDecomposition1() below.
  # However, it has somewhat lower numerical errors and may be used if needed.

  # Adapted from JAMA package, J. Hicklin, C. Moler, P. Webb, R.F. Boisvert,
  # B. Miller, R. Pozo, K. Remington, JAMA: A Java Matrix Package,
  # online at http://math.nist.gov/javanumerics/jama/ (2005),
  # and optimized for Perl.

  my $A = shift;
    # Use a "left-looking", dot-product, Crout/Doolittle algorithm.

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
  my @LUcolj = ();

  # Outer loop.

  for (my $j = 0; $j < $n; $j++) {

    # Make a copy of the j-th column to localize references.

    for (my $i = 0; $i < $n; $i++) {
      $LUcolj[$i] = $LU[$i]->[$j];
    }

    # Apply previous transformations.

    for (my $i = 0; $i < $n; $i++) {
      my $LUrowi = $LU[$i];

      # Most of the time is spent in the following dot product.

      my $kmax = $i < $j ? $i : $j;
      my $s = 0.0;
      for (my $k = 0; $k < $kmax; $k++) {
        $s += $LUrowi->[$k]*$LUcolj[$k];
      }

      $LUrowi->[$j] = $LUcolj[$i] -= $s;
    }

    # Find pivot and exchange if necessary.

    my $p = $j;
    for (my $i = $j+1; $i < $n; $i++) {
      $p = $i if (abs($LUcolj[$i]) > abs($LUcolj[$p]));
    }
    if ($p != $j) {
      my $t = $LU[$p]; $LU[$p] = $LU[$j]; $LU[$j] = $t;
      my $k = $piv[$p]; $piv[$p] = $piv[$j]; $piv[$j] = $k;
      #$pivsign = -$pivsign;
    }

    # Compute multipliers.

    if ( ($j < $n) && $LU[$j]->[$j] ) {
      my $LUjj = 1/$LU[$j]->[$j];
      for (my $i = $j+1; $i < $n; $i++) {
        $LU[$i]->[$j] *= $LUjj;
      }
    } else {
      # Error: matrix is not invertible. Return the index of error row
      return ($j);
    }
  }
  return (\@LU,\@piv);
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
sub getPermutedMatrix($$) {   #12/09/2009 6:44PM
############################################################################
  # Get a submatrix.
  #@param $A   square matrix n x n
  #@param $r   Array of row indices.
  #@return     A(r(:),:)

  my ($A,$r) = @_;
  my $B = [];
  my $num_rows = $#{$r} + 1;
  for (my $i = 0; $i < $num_rows; $i++) {
    $B->[$i] = [];
    push(@{$B->[$i]},@{$A->[$r->[$i]]});
  }
  return $B;
} ##getPermutedMatrix($$)



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

  # Copy right hand side with pivoting
  my $n = $#{$LU} + 1; # number of rows of square matrix $LU and of $B
  my $nx = $#{$B->[0]} + 1; # number of columns of $B
  my $X = getPermutedMatrix($B,$piv);

  # Solve L*Y = B(piv,:)
  for (my $k = 0; $k < $n; $k++) {
    my $Xk = $X->[$k];
    for (my $i = $k+1; $i < $n; $i++) {
      # Take advantage of sparsity of $A: skip the inner cycle
      # if $LU->[$i]->[$k] is zero.
      if ( $LU->[$i]->[$k] ) {
        my $LUik = $LU->[$i]->[$k];
        my $Xi = $X->[$i];
        for (my $j = 0; $j < $nx; $j++) {
          $Xi->[$j] -= $Xk->[$j]*$LUik;
        }
      }
    }
  }
  # Solve U*X = Y;
  for (my $k = $n-1; $k >= 0; $k--) {
    my $Xk = $X->[$k];
    my $LUkk = 1/$LU->[$k]->[$k];
    for (my $j = 0; $j < $nx; $j++) {
      $Xk->[$j] *= $LUkk;
    }
    for (my $i = 0; $i < $k; $i++) {
      # Take advantage of sparsity of $A: skip the inner cycle
      # if $LU->[$i]->[$k] is zero.
      if ( $LU->[$i]->[$k] ) {
        my $LUik = $LU->[$i]->[$k];
        my $Xi = $X->[$i];
        for (my $j = 0; $j < $nx; $j++) {
          $Xi->[$j] -= $Xk->[$j]*$LUik;
        }
      }
    }
  }
  return $X;
} ##LUSolve($$)

############################################################################
sub CholeskyDecomposition($) {
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

      # Transpose the right hand side.
  my $XT = [];
  for ( my $i = 0; $i< $n; $i++ ) {
    for ( my $j = 0; $j< $nx; $j++ ) {
      $XT->[$j]->[$i] = $B->[$i]->[$j];
    }
  }

      # Solve L*Y = X;
  for (my $k = 0; $k < $n; $k++) {
    my $Lk = $L->[$k];
    my $Lkk = $Lk->[$k];
    for (my $j = 0; $j < $nx; $j++) {
      my $XTj = $XT->[$j];
      my $Xkj = $XTj->[$k];
      for (my $i = 0; $i < $k; $i++) {
        $Xkj -= $XTj->[$i]*$Lk->[$i];
      }
      $XTj->[$k] = $Xkj / $Lkk;
    }
  }

    # Solve L'*X = Y;
  for (my $k = $n-1; $k >= 0; $k--) {
    my $Lkk = $L->[$k]->[$k];
    for (my $j = 0; $j < $nx; $j++) {
      my $XTj = $XT->[$j];
      my $Xkj = $XTj->[$k];
      for (my $i = $k+1; $i < $n; $i++) {
        $Xkj -= $XTj->[$i]*$L->[$i]->[$k];
      }
      $XTj->[$k] = $Xkj / $Lkk;
    }
  }

    # Transpose $XT to get the solution
  my $X = [];
  for ( my $i = 0; $i< $n; $i++ ) {
    for ( my $j = 0; $j< $nx; $j++ ) {
      $X->[$i]->[$j] = $XT->[$j]->[$i];
    }
  }
  return $X;
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
  my $dim = $#{$A} + 1;
  my $I = [];
  for ( my $i = 0; $i < $dim; $i++ ) {
    $I->[$i] = [];
    for ( my $j = 0; $j < $dim; $j++ ) {
      $I->[$i]->[$j] = 0;
    }
  }
  for ( my $i = 0; $i < $dim; $i++ ) {
    $I->[$i]->[$i] = 1;
  }

  my $B;

  if ( !$UseCholesky ) {
    my ($LU,$piv) = LUDecomposition1($A);
    if ( ref($LU) ne 'ARRAY' ) {
      # Error: matrix is not invertible. Return the index of the error row
      return ($LU);
    }
    $B = LUSolve($LU,$I,$piv);
  } else {
    my $L = CholeskyDecomposition($A);
    if ( ref($L) ne 'ARRAY' ) {
      # Error: matrix is not invertible. Return the index of the error row
      # or -1 if matrix is asymmetric.
      return ($L);
    }
    $B = SolveCholesky($L,$I);
  }

  my $RMS = 0;
  my $BIG = 0;
  if ($InvChk) {
    for (my $i = 0; $i < $dim; $i++) {
      my $PA = $A->[$i];
      for (my $j = 0; $j < $dim; $j++) {
        my $Aij = ($j == $i) ? -1 : 0;
        for (my $k = 0; $k < $dim; $k++) {
          $Aij += $PA->[$k]*$B->[$k]->[$j];
        }
        $Aij *= $Aij;
        $RMS += $Aij;
        if ($Aij > $BIG) {
          $BIG = $Aij;
        }
      }
    }
    $RMS = sqrt($RMS/(0.5*$dim*($dim+1)));
    $BIG = sqrt($BIG);
  }
  return ($B,$RMS,$BIG);
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

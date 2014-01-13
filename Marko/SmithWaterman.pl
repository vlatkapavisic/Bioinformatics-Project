use Time::HiRes;

$MATCH = 2;
$MISMATCH = -1;
$GAP = -2;


#Max(firstValue, secondValue, ... , lastValue);
#out => maximal value
sub Max{
	my $maxVal = @_[0];
	foreach my $value (@_){
		if($value > $maxVal){
			$maxVal = $value;
		}
	}
	return ($maxVal);
}


#SmithWaterman(firstString, secondString);
#out => indices of the first best matrix element
sub SmithWaterman{
	my $x = @_[0];
	my $y = @_[1];
	my $m = length $x;
	my $n = length $y;
	my @h0 = (0) x ($m+1);
	my @h1 = @h0;
	my $maxValue = 0;
	my $maxI = 0;
	my $maxJ = 0;
	for(my $row = 1; $row < $n+1; $row += 1){
		for(my $col = 1; $col < $m+1; $col += 1){
			$h1[$col] = &Max(0, $h0[$col-1]+ (substr($x, $col-1, 1) eq substr($y, $row-1, 1) ? $MATCH : $MISMATCH), $h1[$col-1] + $GAP, $h0[$col] + $GAP);
			if ($h1[$col] > $maxValue){
				$maxValue = $h1[$col];
				$maxI = $col-1;
				$maxJ = $row-1;
			}
		}
		@h0 = @h1[0..$#h1];
		@h1 = (0) x ($m+1);
	}
	($maxI, $maxJ);
}

#LocalToGlobal(firstString, secondString);
#out => strings that will be used in the global alignment
sub LocalToGlobal{
	my $xLocal = @_[0];
	my $yLocal = @_[1];
	my @border = &SmithWaterman($xLocal, $yLocal);
	my $xReverse = scalar reverse substr ($xLocal, 0, $border[0]+1);
	my $yReverse = scalar reverse substr($yLocal, 0, $border[1]+1);
	@border = &SmithWaterman($xReverse, $yReverse);
	my $xNormal = reverse substr($xReverse, 0, $border[0]+1);
	my $yNormal = reverse substr($yReverse, 0, $border[1]+1);
	
	($xNormal, $yNormal);
}

#NWScore(firstString, secondString);
#out => last row in the NW matrix
sub NWScore{
	my @x = split("", @_[0]);
	my @y = split("", @_[1]);
	my $xLen = @x+1;
	my $yLen = @y+1;
	my @score0 = (0) x ($xLen);
	my @score1 = (0) x ($xLen);
	for($i = 1; $i < $yLen; $i += 1){
		$score0[$i] = $score0[$i-1] + $GAP;
	}
	for($i = 1; $i < $xLen; $i += 1){
		$score1[0] = $score0[0] + $GAP;
		for($j = 1; $j < $yLen; $j += 1){
			$score1[$j] = &Max($score0[$j-1] + ($x[$i-1] eq $y[$j-1] ? $MATCH : $MISMATCH), $score1[$j-1] + $GAP, $score0[$j] + $GAP);
		}
		@score0 = @score1;
		@score1 = (0) x ($xLen);
	}
	return (@score0);
}

#NeedlemanWunsch(firstString, secondString);
#out => optimal alignment using the Needleman-Wunsch algorithm, used only for Nx1 or 1xM
sub NeedlemanWunsch{
	my @x = split("", @_[0]);
	my @y = split("", @_[1]);
	my $xLen = @x;
	my $yLen = @y;
	my $rowLen = $xLen+1;
	my $colLen = $yLen+1;
	my @f = (0) x (($rowLen) * ($colLen)); # f[i][j] = f[i+(xLen+1)*j]
	for($i = 0; $i < $rowLen; $i += 1){
		$f[$i] = $GAP*$i;
	}
	for($i = 0; $i < $colLen; $i += 1){
		$f[$i*($xLen+1)] = $GAP*$i;
	}
	for($i = 1; $i < $rowLen; $i += 1){
		for($j = 1; $j < $colLen; $j += 1){
			$f[$i + $rowLen*$j] = &Max($f[($i-1)+ $rowLen*($j-1)] + ($x[$i-1] eq $y[$j-1] ? $MATCH : $MISMATCH), $f[($i-1)+$rowLen*$j] + $GAP, $f[$i+$rowLen*($j-1)] + $GAP);
		}
	}
	
	my @alx = ();
	my @aly = ();
	
	$i = $xLen;
	$j = $yLen;
	
	while((($i > 0) || ($j > 0))){
		if(($i > 0) && ($j > 0) && ($f[$i+ $rowLen*$j] == $f[$i-1 + $rowLen*($j-1)]+ ($x[$i-1] eq $y[$j-1] ? $MATCH : $MISMATCH))){
			push(@alx, $x[$i-1]);
			push(@aly, $y[$j-1]);
			$i -= 1;
			$j -= 1;
		}
		elsif(($i > 0) && ($f[$i + $rowLen*$j] == $f[$i-1 + $rowLen*$j] + $GAP)){
			push(@alx, $x[$i-1]);
			push(@aly, '-');
			$i -= 1;
		}
		elsif(($j > 0) && ($f[$i + $rowLen*$j] == $f[$i + $rowLen*($j-1)] + $GAP)){
			push(@alx, '-');
			push(@aly, $y[$j-1]);
			$j -= 1;
		}
	}
	
	return (join("", reverse @alx), join("", reverse @aly));
}

#PartitionY(firstRow, secondRow);
#out => index where the sum of elements from both lists is max
sub PartitionY{
	my $sl = @_[0];
	my $sr = @_[1];
	my @sl = @$sl;
	my @sr = @$sr;
	my $lenSr = @sr;
	my @srr = reverse @sr;
	my $biggest = $sl[0] + $srr[0];
	my $index = 0;
	for($i = 1; $i < $lenSr; $i += 1){
		$sum = $sl[$i] + $srr[$i];
		if($sum > $biggest){
			$biggest = $sum;
			$index = $i;
		}
	}
	$index;
}

#Hirschberg(firstString, secondString);
#out => complete optimal global alignment of two strings
sub Hirschberg{
	my $z = "";
	my $w = "";
	my @x = split("", @_[0]);
	my @y = split("", @_[1]);
	my $xLen = @x;
	my $yLen = @y;
	if($xLen == 0 || $yLen == 0){
		if($xLen == 0){
			$z = "";
			$w = "";
			for($i = 0; $i < $yLen; $i += 1){
				$z .= '-';
				$w .= $y[$i];
			}
		}
		elsif($yLen == 0){
			for($i = 0; $i < $xLen; $i += 1){
				$z .= $x[$i];
				$w .= '-';
			}
		}
	}
	elsif($xLen == 1 || $yLen == 1){
		($z, $w) = &NeedlemanWunsch(join("", @x), join("", @y));
	}
	else{
		my $xmid = int($xLen/2);
		my @scoreL = &NWScore(join("", @x[0..$xmid-1]), join("", @y));
		my @scoreR = &NWScore(join("", reverse @x[$xmid..$#x]), join("", reverse @y));
		my $ymid = &PartitionY(\@scoreL, \@scoreR);
		my ($zl, $wl) = &Hirschberg(join("", @x[0..$xmid-1]), join("", @y[0..$ymid-1]));
		my ($zr, $wr) = &Hirschberg(join("", @x[$xmid..$#x]), join("", @y[$ymid..$#y]));
		$z = $zl . $zr;
		$w = $wl . $wr;
	}
	return ($z, $w);
}


#MAIN
{
	#Arguments check
	if(scalar @ARGV < 2){
		print "Correct usage: \"perl $0 <fileName1> <fileName2> <optional:FileOut>\"";
		die "Wrong call";
	}
	if(scalar @ARGV >= 3){
		$fileOut = $ARGV[2];
	}
	else{
		$fileOut = "output.txt";
	}
	#File input
	open FILEONE, @ARGV[0];
	open FILETWO, @ARGV[1];

	while(defined( my $line = <FILEONE>)){
		$line =~ s/^>(.*)$//g;
		$line =~ s/\s//g;
		$stringOne .= $line;
	}

	while(defined( my $line = <FILETWO>)){
		$line =~ s/^>(.*)$//g;
		$line =~ s/\s//g;
		$stringTwo .= $line;
	}

	close FILEONE;
	close FILETWO;
	
	#String error check
	if($stringOne le "" || $stringTwo le ""){
		print "One of the strings is empty!";
		die "Bad input";
	} 

	$startTime = [Time::HiRes::gettimeofday];
	
	#Start of algorithm
	my($globalFirst, $globalSecond) = &LocalToGlobal($stringOne, $stringTwo);
	my($alignedFirst, $alignedSecond) = &Hirschberg($globalFirst, $globalSecond);
	#End of algorithm
	$endTime = [Time::HiRes::gettimeofday];

	$interval = Time::HiRes::tv_interval ($startTime, $endTime);
	
	#File output
	open FILEOUT, ">", $fileOut;
	print FILEOUT ">Hirschberg, sequence lengths ", length $stringOne, " and ", length $stringTwo, ", time taken: $interval s.\n";
	print FILEOUT $alignedFirst, "\n>Second string:\n", $alignedSecond;
	close FILEOUT;
}



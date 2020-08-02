# remisct v1.4, 14 Jan 2015
# Remove duplicate segments in BED (Browser Extensible Display) file
# favoring longest.
# Of two overlapping segments the longer one remained intact and
# the shorter one was trimmed.
# Author: Lev I. Uralsky (Institute of Molecular Genetics, Moscow, Russia)
# v1.0,  22 Dec 2012 - initial release
# v1.1,  22 Dec 2012
# v1.2,   4 Jan 2013
# v1.3,   8 Jan 2013
# v1.3a, 27 Nov 2014 - correct var chrID
# v1.4,  14 Jan 2015 - corrected syntax; added commentaries;
#                      removing unused elements
#
# Usage: gawk -v header=1 -f remisct-v1.3.awk input.bed > output.bed
# Recommended run two times, because second pass may also make changes.

BEGIN {
  # setting input field separator to tab
  FS = "\t";
  # setting output field separator to tab
  OFS = "\t";
  # setting array subscript separator to a semicolon
  SUBSEP = ";";
  # setting default value, from which line to processing begins
  # it may be reset by passing -v chrID=something
#  if(!chrID) chrID = "chr";
  # setting counter default value
  i = 1;
}

header && /^track|^browser/ { print }

{
  # on each line of input file we are
  if ($1 ~ chrID) {
    # getting chromosome id
    chrom = $1;
    # getting the starting position in string 10 character length
    mStart = sprintf("%10s", $2);
    # getting the ending position in string 10 character length
    mEnd = sprintf("%10s", $3);
    # getting name of BED line
    mType = $4;
    # getting score value
    score = $5;
    # getting the strand
    strand = $6;
    # getting field with value, usually, same as starting position
    tStart = $7;
    # getting field with value, usually, same as ending position
    tEnd = $8;
    # getting RGB value, which will determine the display color of the data
    # contained in this line
    itemRGB = $9;
    # setting variable with calculated length of segment
    mLength = mEnd - mStart;
    # get score/length ratio
    scoreToLength = score / mLength;

    # add all received parameters as string with separators to the array,
    # indexed by chromosome id and counter with line number
    aLine[chrom,i] = chrom SUBSEP mStart SUBSEP mEnd SUBSEP mType SUBSEP score SUBSEP strand SUBSEP tStart SUBSEP tEnd SUBSEP itemRGB SUBSEP mLength SUBSEP scoreToLength;
    # second array indexed by chromosome id contains all the relevant line numbers
    aChrs[chrom] = aChrs[chrom] SUBSEP i;

    # removing duplicate coordinates:
    # making new arrays indexed by chr id and starting or ending positions
    # first check the existing values
    if ((chrom SUBSEP mStart) in aStart) {
      # then gets the line from array
      # counter as second index in array we get from the already filtered array
      lineAlreadyIn = aLine[chrom,aStart[chrom,mStart]];
      # using the functions we get the field with length from a delimited string
#      mLengthPrev = GetField(lineAlreadyIn, 10);
      scoreToLengthPrev = GetField(lineAlreadyIn, 11);
      # replacing values in array only if length of current feature greater
      # than exist
#      if (mLength > mLengthPrev) {
      if (scoreToLength > scoreToLengthPrev) {
        aStart[chrom,mStart] = i;
        aEnd[chrom,mEnd] = i;
      }
    }
    # and the same operations for end positions
    if ((chrom SUBSEP mEnd) in aEnd) {
      lineAlreadyIn = aLine[chrom,aEnd[chrom,mEnd]];
#      mLengthPrev = GetField(lineAlreadyIn, 10);
      scoreToLengthPrev = GetField(lineAlreadyIn, 11);
#      if (mLength > mLengthPrev) {
      if (scoreToLength > scoreToLengthPrev) {
        aStart[chrom,mStart] = i;
        aEnd[chrom,mEnd] = i;
      }
    }
    # adds to the array new values, check the non-zero length
    if (!((chrom SUBSEP mStart) in aStart) && !((chrom SUBSEP mEnd) in aEnd) && mLength) {
      aStart[chrom,mStart] = i;
      aEnd[chrom,mEnd] = i;
    }
    # increments the line counter
    i++;
  }
}
# after all lines been procsseed
END {
  # merging arrays with the filtered coordinates
  for (elemE in aEnd) {
    # making a new array, which indexed by line numbers and
    # contains chr id and ending position
    aEndR[aEnd[elemE]] = elemE;
  }
  # now looking to array with starting positions
  for (elemS in aStart) {
    # if both arrays refer to the same line number
    if (aStart[elemS] in aEndR) {
      # getting ending position
      cEnd = aEndR[aStart[elemS]];
      # getting chr id
      chrom = GetField(elemS,1);
      # making a new array indexed by chr id and line numbers,
      # contains starting and ending segment position
      # like this: chr1;10000250;chr1;10001000
      aCoords[chrom,aStart[elemS]] = elemS SUBSEP cEnd;
    }
  }

  # making sorted array in current group (chr id),
  # indexed by the starting position
  # with array which contains line numbers for each chr id
  for (elemChr in aChrs) {
    # removing first semicolon, because first element is empty
    sub(/^;/, "", aChrs[elemChr]);
    # making a new array with line numbers indexed by the order of elements
    split(aChrs[elemChr], aLineInGroup, SUBSEP);
    # and with this array by line numbers in current chr id group
    for (elemLIG in aLineInGroup) {
      # if element with chr id and line number was found in filtered array,
      # which contains starting and ending positions
      if ((elemChr SUBSEP aLineInGroup[elemLIG]) in aCoords) {
        # getting the field with the starting position
        cStart = GetField(aCoords[elemChr,aLineInGroup[elemLIG]], 2);
        # getting the starting position in string 10 character length
        cStart = sprintf("%10s", cStart);
        # making a new array with starting position in current group
        aStartIG[cStart] = aLineInGroup[elemLIG];
      }
    }
    # making a new array sorted by starting position,
    # which indexed by the order of elements
    n = asorti(aStartIG, aStartIGS);
    # and with this array
    for (j=1; j<=n; j++) {
      # getting starting position as number
      currStart = aStartIGS[j]+0;
      # getting the field with ending position
      currEnd = GetField(aCoords[elemChr,aStartIG[aStartIGS[j]]], 4)+0;
      # getting the previous field with ending position
      prevEnd = GetField(aCoords[elemChr,aStartIG[aStartIGS[j-1]]], 4)+0;
      # getting whole line by chr id and starting position
      line = aLine[elemChr,aStartIG[aStartIGS[j]]];
      # removing length placed in last field
#      gsub(/[[:digit:]]+$/, "", line);
      sub(/[[:digit:]]+.[[:digit:]\.]+$/, "", line);
      # removing extra spaces
      gsub(/[[:space:]]/, "", line);
      # searching overlapping elements
      if (currStart > prevEnd) {
        # replacing subscript separators by output field separator
        gsub(SUBSEP, OFS, line);
        # replacing extra spaces at the end of line
        gsub(/[[:space:]]$/, "", line);
        # printing the resulting line
        print line;
      } else {
        if (currEnd > prevEnd) {
          # changing starting position in the next line
          gsub(SUBSEP currStart SUBSEP, SUBSEP prevEnd+1 SUBSEP, line);
          # replacing subscript separators by output field separator 
          gsub(SUBSEP, OFS, line);
          # replacing extra spaces at the end of line
          gsub(/[[:space:]]$/, "", line);
          # printing the resulting line
          print line;
        }
      }
    }
    # deleting array for current chr id group
    delete aStartIG;
  }
}

# function which helps get the field from a delimited string
function GetField(string,field,  n) {
  if(length(string"") > 0) {
    n = split(string, array, SUBSEP);
    if (field > 0) {
      if (field <= n) {
        result = array[field];
      }
    } else {
      result = string;
    }
  } else {
    return string;
  }
  return result;
}

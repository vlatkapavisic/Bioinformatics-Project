import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class Program {

	/**
	 * @param args
	 */
	public static int gapPen = -1;
	public static int matchPen = +2;
	public static int mismatchPen = -1;
	
	public static void main(String[] args) throws IOException {
		startWithReading(args);
	}

	// Finds the alignment of two strings which it'll read from a file
	public static void startWithReading(String[] args) throws IOException {
		StringPair ab = readFromFile(args); 			// Gets the input strings from the given files
		
		long startTime = System.nanoTime(); 					// Starts measuring time
	    ab = SmithWaterman(ab.a, ab.b); 						// Finds the local alignment of the input strngs
		ab = Hirschberg(ab.a, ab.b); 							// Starts the recursion which finds the global alignment of the strings
		long endTime = System.nanoTime();						// Finishes measuring time
		System.out.println((endTime - startTime) + "ns");		// Displays the time (in nanosecond)s it took to find the alignment
		
		printResult(ab.a, ab.b, args);							// Displays the result either in the console or writes it in a file
	}

	// Reads the input strings from two separate files
	public static StringPair readFromFile(String[] args) throws IOException {
		String aPath, bPath;
		if(args.length < 2) {			// If the file paths were not given as arguments, it asks the user for the location of the files
			System.out.println("Enter the path of the file containing the first sequence:");
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
	    	aPath = br.readLine();

			System.out.println("Enter the path of the file containing the second sequence:");
	    	bPath = br.readLine();
		}
		else {
			aPath = args[0];
			bPath = args[1];
		}
	    
	    String a = readFile(aPath); 	// Reads the input files
	    String b = readFile(bPath); 
	    return new StringPair(a, b);
	}
	
	// Displays the final result and writes it to a file if the user wants to
	private static void printResult(String a, String b, String[] args) throws IOException {
		String mid = "";
		for(int i = 0; i < a.length(); i++) {
			if (a.charAt(i) == b.charAt(i)) mid += "|";							// An easier way to display if the characters match or
			else if(a.charAt(i) == '-' || b.charAt(i) == '-') mid += " ";		// if there's an empty space in one of the alignments
			else mid += "x";
		}
		
		String outPath;
		if(args.length < 3) {					// If the output path is not given as an arugment, the program asks the user to write it in the console
			System.out.println("Enter the path of the output file or press ENTER if you don't want to save it to a file:");
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			outPath = br.readLine();
		}
		else outPath = args[2];
		
	    if (!outPath.equals("")) {				// If the output file is set, we write the alignment in it
	        BufferedWriter out = new BufferedWriter(new FileWriter(outPath));
	        out.write(">Optimal local alignment of the first sequence:\r\n");
	        out.write(a + "\r\n");
	        out.write(">Optimal local alignment of the second sequence:\r\n");
	        out.write(b + "\r\n");
	        out.close();
	    }
	    else {									// If the output file is not set, the alignment is diplayed in the console
			System.out.println(a);
			System.out.println(mid);
			System.out.println(b);
	    }
	}

	// Uses Hirschberg's algorithm to recursively find the alignment of the two strings
	private static StringPair Hirschberg(String a, String b) {
		String z = "";
		String w = "";
		if (a.length() == 0 || b.length() == 0) {				// If one of the strings is empty, we align the other string's characters
	    	if (a.length() == 0) {								// with empty spaces in that one
	    		z = "";
	    		w = "";
	    		for (int i = 0; i < b.length(); i++) {
	            	z += '-';
	               	w += b.charAt(i);
	    		}
	    	}
	    	else if (b.length() == 0) {
	    		for (int i = 0; i < a.length(); i++) {
	    			z += a.charAt(i);
	    			w += '-';
	    		}
	    	}
		}
	    else if (a.length() == 1 || b.length() == 1) {			// If one of the strings has a length of 1, we can use the normal Needleman-Wunsch
	    	StringPair zw = NeedlemanWunsch(a, b);				// algorithm to ailgn them without the memory requirements exceeding linear
	    	z = zw.a; w = zw.b;
	    }
	    else {
	    	int xmid = a.length()/2; 														// Otherwise we find the middle x-value and use it to
	        int[] score_l = NWScore(a.substring(0, xmid), b);								// find the middle y-value through the Needleman-Wunsch 
	        int[] score_r = NWScore(reverseString(a.substring(xmid)), reverseString(b));	// algorithm by finding the maximum value of the middle edge
	        int ymid = PartitionY(score_l, score_r);										// between the middle two rows and then we align the
	        StringPair zlwl = Hirschberg(a.substring(0, xmid), b.substring(0, ymid));		// top-left and bottom-right subsequences recursively
	        StringPair zrwr = Hirschberg(a.substring(xmid), b.substring(ymid));
	        z = zlwl.a + zrwr.a;
	        w = zlwl.b + zrwr.b;
	    }
		return new StringPair(z, w);
	}
	
	// Finds the last row using the Needleman-Wunsch algorithm
	private static int[] NWScore(String a, String b) {
		ArrayList<Integer> row0 = new ArrayList<Integer>();
		ArrayList<Integer> row1 = new ArrayList<Integer>();
		row0.add(0);												// Creates the first row of the table by giving 
		for (int j = 1; j < b.length() + 1; j++) {					// each element the value of index*gapPenalty
			row0.add(row0.get(j - 1) + gapPen);
		}
		for (int i = 1; i < a.length() + 1; i++) {					// Uses the row0 to find the row1 (bottom row) by finding the best way to align
			row1.add(row0.get(0) + gapPen);							// the sequences
			for(int j = 1; j < b.length() + 1; j++) {
				row1.add(max(row0.get(j - 1) + matchPenalty(a.charAt(i - 1), b.charAt(j - 1)), row1.get(j - 1) + gapPen, row0.get(j) + gapPen));
			}
			row0 = row1;											// Makes the former bottom row the top row and creates a new bottom row for the next
			row1 = new ArrayList<Integer>();						// step o the for loop
		}
		return listToArray(row0);
	}

	// Creates an integer array from an integer list
	private static int[] listToArray(ArrayList<Integer> l) {
		int[] ret = new int[l.size()];
		for (int i = 0; i < l.size(); i++) {
			ret[i] = l.get(i);
		}
		return ret;
	}

	// Finds the index where the sum of the elements in the two rows is the biggest
	private static int PartitionY(int[] sl, int[] sr) {
		ArrayList<Integer> sum_list = new ArrayList<Integer>();
		int[] srr = reverseArray(sr);
		for (int i = 0; i < sl.length; i++) {
			sum_list.add(sl[i] + srr[i]);
		}
		return maxList(sum_list);
	}

	// Finds the element with the highest value in a list of integers and returns it's index
	private static int maxList(ArrayList<Integer> l) {
		int max = l.get(0), maxI = 0;
		for(int i = 1; i < l.size(); i++) {
			if(l.get(i) > max) {
				max = l.get(i);
				maxI = i;
			}
		}
		return maxI;
	}

	// Uses the Needleman-Wunsch algorithm to find the alignment when at least one of
	// the strings has the size of 1 and therefore it uses linear space requirements
	private static StringPair NeedlemanWunsch(String a, String b) {
		int[][] f = new int[a.length() + 1][b.length() + 1];
		for(int i = 0; i < a.length() + 1; i++) {
			f[i][0] = gapPen*i;
		}
		for(int j = 0; j < b.length() + 1; j++) {
			f[0][j] = gapPen*j;
		}
		for(int i = 1; i < a.length() + 1; i++) {
			for(int j = 1; j < b.length() + 1; j++) {
				f[i][j] = max(f[i - 1][j - 1] + matchPenalty(a.charAt(i - 1), b.charAt(j - 1)), f[i - 1][j] + gapPen, f[i][j - 1] + gapPen);
			}
		}
		String alx = "";
		String aly = "";
		int i = a.length();
		int j = b.length();
		
		while(i > 0 || j > 0) {
			if(i > 0 && j > 0 && f[i][j] == f[i - 1][j - 1] + matchPenalty(a.charAt(i - 1), b.charAt(j - 1))) {
				alx += a.charAt(i - 1);
				aly += b.charAt(j - 1);
				i--; j--;
			}
			else if (i > 0 && f[i][j] == f[i - 1][j] + gapPen) {
				alx += a.charAt(i - 1);
				aly += '-';
				i--;
			}
			else if (j > 0 && f[i][j] == f[i][j - 1] + gapPen) {
				alx += '-';
				aly += b.charAt(j - 1);
				j--;
			}
		}
		return new StringPair(reverseString(aly), reverseString(alx));
	}

	// Finds the local alignment of the two strings using the Smith-Waterman algorithm
	private static StringPair SmithWaterman(String a, String b) {
		StringPair temp = reduceStrings(a, b);
        a = reverseString(temp.a);
        b = reverseString(temp.b);

        temp = reduceStrings(a, b);
        a = reverseString(temp.a);
        b = reverseString(temp.b);
        
        if(a.length() < b.length()) {
        	String temp2 = a;
        	a = b;
        	b = temp2;
        }
        
        return new StringPair(a, b);
	}
	
	// Reverses the given string
	private static String reverseString(String a) {
        return new StringBuilder(a).reverse().toString();
    }
	
	// Uses the Smith-Waterman algorithm to reduce the length of the strings depending on their global alignment
	private static StringPair reduceStrings(String a, String b) {
		int[] row1 = new int[a.length() + 1];
        int[] row2 = new int[a.length() + 1];
        int maxValue = 0, maxI = 0, maxJ = 0;
        int left, up, diag;

        for (int i = 0; i < b.length(); i++)			// Goes through every element in the table in order to find the one with the maximum value
        {
            for (int j = 1; j < row1.length; j++)
            {
                left = row2[j - 1] + gapPen; if (left < 0) left = 0;
                up = row1[j] + gapPen; if (up < 0) up = 0;
                diag = row1[j - 1] + matchPenalty(a.charAt(j - 1), b.charAt(i)); if (diag < 0) diag = 0;

                row2[j] = max(left, up, diag);
                if (row2[j] > maxValue)
                {
                    maxValue = row2[j];
                    maxI = i;
                    maxJ = j;
                }
            }

            row1 = row2;
            row2 = new int[a.length() + 1];
        }
        
        return new StringPair(a.substring(0, maxJ), b.substring(0, maxI + 1));	// And in the end it cuts the strings at the index of the max value
	}
	
	// Finds the biggest of the three integers
	private static int max(int left, int up, int diag) {
        if (diag >= left && diag >= up) return diag;
        else if (up >= left) return up;
        else return left;
    }
	
	// Returns the match or mismatch penalty depending on whether the characters are equal
	private static int matchPenalty(char a, char b) {
        return (a == b) ? matchPen : mismatchPen;
    }
	
	// Reverses the given array
	private static int[] reverseArray(int[] a) {
		int[] ret = new int[a.length];
		for(int i = 0; i < a.length; i++) {
			ret[i] = a[a.length - 1 - i];
		}
		return ret;
	}
	
	// Reads the whole given file and returns it as a string, leaves out the lines starting with ">"
		private static String readFile(String path) throws IOException {
			String content = "";
			BufferedReader br = new BufferedReader(new FileReader(path));
			String line;
			while ((line = br.readLine()) != null) {
			   if(line.charAt(0) != '>') content += line;
			}
			br.close();
			return content;
		}
		
	// A class that's used to hold two strings, so it's easier to return them from a function
	public static class StringPair {
		public String a;
		public String b;
		
		public StringPair (String _a, String _b) {
			a = _a;
			b = _b;
		}
	}
}

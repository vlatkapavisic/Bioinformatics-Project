using System;
using System.IO;
using System.Text;

namespace Hirschberg
{
    class Program
    {
        private static void Main(string[] args)
        {
            try
            {
                //TimeTesting.TestHirschberg();
                //Console.ReadLine();
                //return;
                FileParser.ReadInputs(args);
                char[] x = FileParser.Sequence1.ToCharArray(); 
                //char[] x = Console.ReadLine().ToCharArray();
                char[] y = FileParser.Sequence2.ToCharArray();
                //char[] y = Console.ReadLine().ToCharArray();
               
                bool isWritingToConsole = false;
                bool isCheckingMemory = false;
                bool isCheckingTime = false;

                for (int i = 0; i < args.Length; i++)
                {
                    if (args[i] == "-c")
                    {
                        isWritingToConsole = true;
                    }
                    else if (args[i] == "-m")
                    {
                        isCheckingMemory = true;
                    }
                    else if(args[i] == "-t")
                    {
                        isCheckingTime = true;
                    }
                }
                
                DateTime startTime = DateTime.Now;
                
                SequenceAlignment.Align(x, y, isWritingToConsole, isCheckingMemory);
                if(isCheckingTime)
                    Console.WriteLine("Time taken: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");
                if(isCheckingMemory)
                    Console.WriteLine("Memory taken: " + SequenceAlignment.MemoryUsageMax + "bytes");
                FileParser.TryStoringToFileLessMemory();
                Console.WriteLine("Finished successfully");
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }

            //Console.ReadLine();
        }
		
            /// <summary>
        /// First sequence
        /// </summary>
        public static string Sequence1 { get; private set; }

        /// <summary>
        /// Second sequence
        /// </summary>
        public static string Sequence2 { get; private set; }

        private static string OutputFilePath;

        /// <summary>
        /// Input arguments should be one of the following:
        /// 1. version
        /// 2 relative or absolute paths to FASTA or regular txt files containing one input sequence each
        /// example: programName absolutePath1 absolutePath2
        ///          programName relativePath1 relativePath2
        /// file extension is irrelevant
        /// 
        /// 2. version
        /// 1 path containing both input sequences
        /// example: programName absolutePath1
        ///          programName relativePath1
        /// 
        /// files are assumed to be FASTA files if first line starts with '>'
        /// </summary>
        /// <param name="args">Program input arguments</param>
        public static void ReadInputs(string[] args)
        {
            // if memory error occurs or there are no arguments, input sequences cannot be read
            if (args == null || args.Length < 2)
                throw new ArgumentException("Invalid input arguments, see README.txt for help.");
            
            string filePath1 = CheckFileExists(args[0]);
            string filePath2 = CheckFileExists(args[1]);
            ReadInputsFrom2Files(filePath1, filePath2);

            if (args.Length == 3)
            {
                string directoryPath = Path.GetDirectoryName(args[2]);
                if (!Directory.Exists(directoryPath))
                {
                    directoryPath = AppDomain.CurrentDomain.BaseDirectory + args[2];
                    if (!Directory.Exists(Path.GetDirectoryName(directoryPath)))
                        throw new ArgumentException(
                            "Third argument for output path is invalid, see README.txt for help.");
                }
                OutputFilePath = directoryPath;
            }
        }

        /// <summary>
        /// Checks relative and apsolute path
        /// </summary>
        /// <param name="path">input path for testing</param>
        /// <returns>string that represents path where file exists, or throws FileNotFoundException</returns>
        private static string CheckFileExists(string path)
        {
            string relativePathPrefix = AppDomain.CurrentDomain.BaseDirectory;

            if (!File.Exists(relativePathPrefix + path))
            {
                if (!File.Exists(path))
                    // if any of given files does not exist, input sequences cannot be aligned
                    throw new FileNotFoundException("File not found. Searched paths: " + path + "; " +
                                                    relativePathPrefix + path);
                else
                    return path;
            }
            else 
                return relativePathPrefix + path;
        }

        /// <summary>
        /// Reads files from 2 input paths, checked that file exists
        /// </summary>
        /// <param name="path1">path for the first file</param>
        /// <param name="path2">path for the second file</param>
        private static void ReadInputsFrom2Files(string path1, string path2)
        {
            StringBuilder builder = new StringBuilder();
            StreamReader reader = null;
            try
            {
                reader = new StreamReader(path1);

                string line;

                while (!reader.EndOfStream)
                {
                    line = reader.ReadLine().Trim();
                    if (line != null && line.Length != 0 && line[0] != '>')
                        builder.Append(CleanString(line));
                }
                Sequence1 = builder.ToString();

                reader.Close();
                reader = new StreamReader(path2);
                builder.Clear();

                while (!reader.EndOfStream)
                {
                    line = reader.ReadLine().Trim();
                    if(line!=null && line.Length != 0 && line[0] != '>')
                        builder.Append(CleanString(line));
                }
                Sequence2 = builder.ToString();

                if(Sequence1==null || Sequence1.Length == 0 || Sequence2==null || Sequence2.Length == 0)
                    throw new ArgumentException("Error while reading sequences");
            }
            catch (Exception e)
            {
                throw e;
            }
            finally
            {
                if(reader != null)
                    reader.Close();
            }
        }
        
        /// <summary>
        /// Remove invalid ASCII characters from line and return it
        ///  </summary>
        /// <param name="line">Input line from file</param>
        /// <returns>"Clean" string, throws exception if result string is empty</returns>
        private static string CleanString(string line)
        {
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < line.Length; i++)
            {
                if (line[i] > ((char) 32) && line[i] < 127)
                    builder.Append(line[i]);
            }

            line = builder.ToString();

            if(line.Length==0)
                throw new ArgumentException("Line " + line + " did not contain valid characters");

            return line;
        }
        
        /// <summary>
        /// Method that stores alignment data to output file
        /// </summary>
        public static void TryStoringToFile()
        {
            if (string.IsNullOrEmpty(OutputFilePath))
                return;

            if (SequenceAlignment.XResult == null || SequenceAlignment.Matching == null ||
                SequenceAlignment.YResult == null)
                return;

            StreamWriter writer = new StreamWriter(OutputFilePath);
            writer.WriteLine(SequenceAlignment.Score);
            writer.WriteLine(SequenceAlignment.XResult);
            writer.WriteLine(SequenceAlignment.Matching);
            writer.WriteLine(SequenceAlignment.YResult);
            writer.Close();
        }

        public static void TryStoringToFileLessMemory()
        {
            if (string.IsNullOrEmpty(OutputFilePath))
                return;

            if (SequenceAlignment.XResult == null || SequenceAlignment.Matching == null ||
                SequenceAlignment.YResult == null)
                return;

            StreamWriter writer = new StreamWriter(OutputFilePath);
            writer.WriteLine(SequenceAlignment.XResult);
            
            int score = 0;

            char[] x = SequenceAlignment.XResult;
            char[] y = SequenceAlignment.YResult;

            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] == '-' || y[i] == '-')
                {
                    score += SequenceAlignment.GAP;
                    writer.Write(' ');
                }
                else
                {
                    if (x[i] == y[i])
                    {
                        score += SequenceAlignment.MATCH;
                        writer.Write('|');
                    }
                    else
                    {
                        score += SequenceAlignment.MISMATCH;
                         writer.Write('x');
                    }
                }
            }
            writer.WriteLine();
            writer.WriteLine(score);
            writer.WriteLine(SequenceAlignment.YResult);
            writer.Close();
        }
    }
	

	/// <summary>
    /// Class for reading sequences from input files
    /// </summary>
    static class FileParser
    {
        /// <summary>
        /// First sequence
        /// </summary>
        public static string Sequence1 { get; private set; }

        /// <summary>
        /// Second sequence
        /// </summary>
        public static string Sequence2 { get; private set; }

        private static string OutputFilePath;

        /// <summary>
        /// Input arguments should be one of the following:
        /// 1. version
        /// 2 relative or absolute paths to FASTA or regular txt files containing one input sequence each
        /// example: programName absolutePath1 absolutePath2
        ///          programName relativePath1 relativePath2
        /// file extension is irrelevant
        /// 
        /// 2. version
        /// 1 path containing both input sequences
        /// example: programName absolutePath1
        ///          programName relativePath1
        /// 
        /// files are assumed to be FASTA files if first line starts with '>'
        /// </summary>
        /// <param name="args">Program input arguments</param>
        public static void ReadInputs(string[] args)
        {
            // if memory error occurs or there are no arguments, input sequences cannot be read
            if (args == null || args.Length < 2)
                throw new ArgumentException("Invalid input arguments, see README.txt for help.");
            
            string filePath1 = CheckFileExists(args[0]);
            string filePath2 = CheckFileExists(args[1]);
            ReadInputsFrom2Files(filePath1, filePath2);

            if (args.Length == 3)
            {
                string directoryPath = Path.GetDirectoryName(args[2]);
                if (!Directory.Exists(directoryPath))
                {
                    directoryPath = AppDomain.CurrentDomain.BaseDirectory + args[2];
                    if (!Directory.Exists(Path.GetDirectoryName(directoryPath)))
                        throw new ArgumentException(
                            "Third argument for output path is invalid, see README.txt for help.");
                }
                OutputFilePath = directoryPath;
            }
        }

        /// <summary>
        /// Checks relative and apsolute path
        /// </summary>
        /// <param name="path">input path for testing</param>
        /// <returns>string that represents path where file exists, or throws FileNotFoundException</returns>
        private static string CheckFileExists(string path)
        {
            string relativePathPrefix = AppDomain.CurrentDomain.BaseDirectory;

            if (!File.Exists(relativePathPrefix + path))
            {
                if (!File.Exists(path))
                    // if any of given files does not exist, input sequences cannot be aligned
                    throw new FileNotFoundException("File not found. Searched paths: " + path + "; " +
                                                    relativePathPrefix + path);
                else
                    return path;
            }
            else 
                return relativePathPrefix + path;
        }

        /// <summary>
        /// Reads files from 2 input paths, checked that file exists
        /// </summary>
        /// <param name="path1">path for the first file</param>
        /// <param name="path2">path for the second file</param>
        private static void ReadInputsFrom2Files(string path1, string path2)
        {
            StringBuilder builder = new StringBuilder();
            StreamReader reader = null;
            try
            {
                reader = new StreamReader(path1);

                string line;

                while (!reader.EndOfStream)
                {
                    line = reader.ReadLine().Trim();
                    if (line != null && line.Length != 0 && line[0] != '>')
                        builder.Append(CleanString(line));
                }
                Sequence1 = builder.ToString();

                reader.Close();
                reader = new StreamReader(path2);
                builder.Clear();

                while (!reader.EndOfStream)
                {
                    line = reader.ReadLine().Trim();
                    if(line!=null && line.Length != 0 && line[0] != '>')
                        builder.Append(CleanString(line));
                }
                Sequence2 = builder.ToString();

                if(Sequence1==null || Sequence1.Length == 0 || Sequence2==null || Sequence2.Length == 0)
                    throw new ArgumentException("Error while reading sequences");
            }
            catch (Exception e)
            {
                throw e;
            }
            finally
            {
                if(reader != null)
                    reader.Close();
            }
        }
        
        /// <summary>
        /// Remove invalid ASCII characters from line and return it
        ///  </summary>
        /// <param name="line">Input line from file</param>
        /// <returns>"Clean" string, throws exception if result string is empty</returns>
        private static string CleanString(string line)
        {
            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < line.Length; i++)
            {
                if (line[i] > ((char) 32) && line[i] < 127)
                    builder.Append(line[i]);
            }

            line = builder.ToString();

            if(line.Length==0)
                throw new ArgumentException("Line " + line + " did not contain valid characters");

            return line;
        }
        
        /// <summary>
        /// Method that stores alignment data to output file
        /// </summary>
        public static void TryStoringToFile()
        {
            if (string.IsNullOrEmpty(OutputFilePath))
                return;

            if (SequenceAlignment.XResult == null || SequenceAlignment.Matching == null ||
                SequenceAlignment.YResult == null)
                return;

            StreamWriter writer = new StreamWriter(OutputFilePath);
            writer.WriteLine(SequenceAlignment.Score);
            writer.WriteLine(SequenceAlignment.XResult);
            writer.WriteLine(SequenceAlignment.Matching);
            writer.WriteLine(SequenceAlignment.YResult);
            writer.Close();
        }

        public static void TryStoringToFileLessMemory()
        {
            if (string.IsNullOrEmpty(OutputFilePath))
                return;

            if (SequenceAlignment.XResult == null || SequenceAlignment.Matching == null ||
                SequenceAlignment.YResult == null)
                return;

            StreamWriter writer = new StreamWriter(OutputFilePath);
            writer.WriteLine(SequenceAlignment.XResult);
            
            int score = 0;

            char[] x = SequenceAlignment.XResult;
            char[] y = SequenceAlignment.YResult;

            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] == '-' || y[i] == '-')
                {
                    score += SequenceAlignment.GAP;
                    writer.Write(' ');
                }
                else
                {
                    if (x[i] == y[i])
                    {
                        score += SequenceAlignment.MATCH;
                        writer.Write('|');
                    }
                    else
                    {
                        score += SequenceAlignment.MISMATCH;
                         writer.Write('x');
                    }
                }
            }
            writer.WriteLine();
            writer.WriteLine(score);
            writer.WriteLine(SequenceAlignment.YResult);
            writer.Close();
        }
    }
	
	 /// <summary>
    /// Class that implements local alignment with Smith-Waterman cropping and SequenceAlignment using Needleman-Wunsch 
    /// </summary>
    public static class SequenceAlignment
    {
        /// <summary>
        /// Variable that stores aligned result of first sequence
        /// </summary>
        public static char[] XResult;

        /// <summary>
        /// Variable that stores aligned result of second sequence
        /// </summary>
        public static char[] YResult;

        /// <summary>
        /// Variable that stores matching result
        /// </summary>
        public static char[] Matching;

        /// <summary>
        /// Variable that stores alignment score
        /// </summary>
        public static int Score;

        /// <summary>
        /// The price reward for aligning two same characters
        /// </summary>
        public const int MATCH = 2;

        /// <summary>
        /// The price punishment for aligning two different characters
        /// </summary>
        public const int MISMATCH = -1;

        /// <summary>
        /// The price for aligning a char with vertical or horizontal alignment
        /// </summary>
        public const int GAP = -2;

        /// <summary>
        /// Returns the cost of matching given chars
        /// </summary>
        /// <param name="a">char from first sequence</param>
        /// <param name="b">char from second sequence</param>
        /// <returns>The matching price</returns>
        private static int Match(char a, char b)
        {
            return (a == b) ? MATCH : MISMATCH;
        }

        /// <summary>
        /// Calculates SmithWaterman matrix using only 2 rows of memory at a time
        /// and returns the index for x and y coordinate of maximum score
        /// </summary>
        /// <param name="x">First sequence</param>
        /// <param name="y">Second sequence</param>
        /// <param name="max_i">Index for maximum score position on x axis</param>
        /// <param name="max_j">Index for maximum score position on y axis</param>
        private static void SmithWaterman(char[] x, char[] y, out int max_i, out int max_j, out int max_value)
        {
            int m = x.Length;
            int n = y.Length;
            int[] h0 = new int[m + 1];
            int[] h1 = new int[m + 1];
            max_value = 0;
            max_i = 0;
            max_j = 0;
            for (int j = 1; j < n + 1; j++) // column index
            {
                for (int i = 1; i < m + 1; i++) // row index
                {
                    int max = 0;
                    int a = h0[i - 1] + ((x[i - 1] == y[j - 1]) ? MATCH : MISMATCH);
                    if (a > max) max = a;
                    a = h1[i - 1] + GAP;
                    if (a > max) max = a;
                    a = h0[i] + GAP;
                    if (a > max) max = a;

                    h1[i] = max;

                    //h1[i] = Math.Max(Math.Max(0, h0[i - 1] + ((x[i - 1] == y[j - 1]) ? MATCH : MISMATCH)),
                    //Match(x[i - 1], y[j - 1])),
                    //                 Math.Max(h1[i - 1] + GAP, h0[i] + GAP));
                    
                    if (h1[i] > max_value)
                    {
                        max_value = h1[i];
                        max_i = i - 1;
                        max_j = j - 1;
                    }
                }
                h0 = h1;
                h1 = new int[m + 1];
            }
        }
        
        /// <summary>
        /// Subsamples source sequences with borders on optimal 
        /// score positions in SmithWaterman alignment acquired from start and end
        /// </summary>
        /// <param name="x">First sequence</param>
        /// <param name="y">Second sequence</param>
        private static void LocalToGlobal(char[] x, char[] y, out char[] xLocal, out char[] yLocal, out int maxScore)
        {
            int border0, border1;
            SmithWaterman(x, y, out border0, out border1, out maxScore);
            char[] xx = new char[border0 + 1];
            Array.Copy(x, xx, border0 + 1);
            Array.Reverse(xx);
            char[] yy = new char[border1 + 1];
            Array.Copy(y, yy, border1 + 1);
            Array.Reverse(yy);
            SmithWaterman(xx, yy, out border0, out border1, out maxScore);
            xLocal = new char[border0 + 1];
            Array.Copy(xx, xLocal, border0 + 1);
            Array.Reverse(xLocal);
            yLocal = new char[border1 + 1];
            Array.Copy(yy, yLocal, border1 + 1);
            Array.Reverse(yLocal);
        }

        /// <summary>
        /// Returns the last row of Needleman-Wunsch table
        /// </summary>
        /// <param name="x">First sequence</param>
        /// <param name="y">Second sequence</param>
        /// <returns>The last row</returns>
        private static int[] NWScore(char[] x, char[] y)
        {
            int[] score0 = new int[y.Length + 1];
            for (int i = 1; i < y.Length + 1; i++)
                score0[i] = score0[i - 1] + GAP;
            int[] score1 = new int[y.Length + 1];

            for (int i = 1; i < x.Length + 1; i++)
            {
                score1[0] = score0[0] + GAP;
                for (int j = 1; j < y.Length + 1; j++)
                {
                    score1[j] = Math.Max(score0[j - 1] + 
                        ((x[i - 1] == y[j - 1]) ? MATCH : MISMATCH),
                        //Match(x[i - 1], y[j - 1]),
                                        Math.Max(score1[j - 1] + GAP, score0[j] + GAP));
                }
                score0 = score1;
                score1 = new int[y.Length + 1];
            }
            return score0;
        }

        /// <summary>
        /// Calculates the alignment for the simplest cases of Hirschberg
        /// </summary>
        /// <param name="x">The first sequence</param>
        /// <param name="y">The second sequence</param>
        /// <param name="seq2">The first output aligned sequence, containing '-' as gaps</param>
        /// <param name="seq1">The seconde output aligned sequence, containing '-' as gaps</param>
        private static void NeedlemanWunsch(char[] x, char[] y, out char[] seq2, out char[] seq1)
        {
            int i, j;
            int[,] f = new int[x.Length+1,y.Length+1];
            for (i = 0; i < x.Length + 1; i++)
                f[i, 0] = GAP*i;
            for (j = 0; j < y.Length + 1; j++)
                f[0, j] = GAP*j;
            for(i=1; i<x.Length + 1; i++)
                for (j = 1; j < y.Length + 1; j++)
                {
                    int max = f[i - 1, j - 1] + ((x[i - 1] == y[j - 1]) ? MATCH : MISMATCH);
                    int a = f[i - 1, j] + GAP;
                    if (a > max) max = a;
                    a = f[i, j - 1] + GAP;
                    if (a > max) max = a;

                    f[i, j] = max;
                    //Math.Max(f[i - 1, j - 1] + ((x[i - 1] == y[j - 1]) ? MATCH : MISMATCH),
                        //Match(x[i - 1], y[j - 1]),
                      //                  Math.Max(f[i-1,j] + GAP, f[i,j-1] + GAP));
                }

            StringBuilder alx = new StringBuilder();
            StringBuilder aly = new StringBuilder();

            i = x.Length;
            j = y.Length;

            while (i > 0 || j > 0)
            {
                if (i > 0 && j > 0 && f[i, j] ==
                    (f[i - 1, j - 1] + ((x[i - 1] == y[j - 1]) ? MATCH : MISMATCH)))
                    //Match(x[i - 1], y[j - 1])))
                {
                    alx.Append(x[i - 1]);
                    aly.Append(y[j - 1]);
                    i--;
                    j--;
                }
                else if(i > 0 && f[i,j] == f[i-1,j] + GAP)
                {
                    alx.Append(x[i - 1]);
                    aly.Append("-");
                    i--;
                }
                else if (j > 0 && f[i, j] == f[i, j - 1] + GAP)
                {
                    alx.Append("-");
                    aly.Append(y[j - 1]);
                    j--;
                }
            }

            seq1 = alx.ToString().ToCharArray();
            Array.Reverse(seq1);
            seq2 = aly.ToString().ToCharArray();
            Array.Reverse(seq2);
        }
    
        /// <summary>
        /// Calculation of optimal position to cross between two score rows
        /// </summary>
        /// <param name="sl">The first score row, calculated from top</param>
        /// <param name="sr">The second score row, calculated from reverse strings - from the bottom</param>
        /// <returns>The index of maximum score cell</returns>
        private static int PartitionY(int[] sl, int[] sr)
        {
            int max_value = int.MinValue;
            int index = -1;
            int[] srr = new int[sr.Length];
            Array.Copy(sr,srr,sr.Length);
            Array.Reverse(srr);
            for (int i = 0; i < sl.Length; i++)
            {
                int value = sl[i] + srr[i];
                if (value > max_value)
                {
                    max_value = value;
                    index = i;
                }
            }
            return index;
        }

        /// <summary>
        /// Calculates alignment by Hirschberg algorithm using Needleman-Wunsch for finding the middle edge
        /// </summary>
        /// <param name="x">The first sequence</param>
        /// <param name="y">The second sequence</param>
        /// <param name="z">The first sequence result alignment</param>
        /// <param name="w">The second sequence result alignment</param>
        private static void Hirschberg(char[] x, char[] y, out char[] z, out char[] w)
        {
            if (x.Length == 0 || y.Length == 0)
            {
                //StringBuilder zBuilder = new StringBuilder();
                //StringBuilder wBuilder = new StringBuilder();

                if (x.Length == 0)
                {
                    z = new char[y.Length];
                    w = new char[y.Length];
                    Array.Copy(y,w,y.Length);

                    for (int i = 0; i < y.Length; i++)
                    {
                        //zBuilder.Append("-");
                        z[i] = '-';
                        //wBuilder.Append(y[i]);
                    }
                }
                else //if(y.Length == 0)
                {
                    z = new char[x.Length];
                    w = new char[x.Length];
                    Array.Copy(x, z, x.Length);

                    for (int i = 0; i < x.Length; i++)
                    {
                        //zBuilder.Append(x[i]);
                        //wBuilder.Append("-");
                        w[i] = '-';
                    }
                }

                //z = zBuilder.ToString().ToCharArray();
                //w = wBuilder.ToString().ToCharArray();
            }
            else if (x.Length == 1 || y.Length == 1)
            {
                NeedlemanWunsch(x,y,out z, out w);
            }
            else
            {
                char[] zl, wl, zr, wr;

                int xmid = x.Length/2;
                int[] score_l = NWScore(SubCharArray(x,0,xmid,false), y);
                int[] score_r = NWScore(SubCharArray(x, xmid, x.Length, true), SubCharArray(y, 0, y.Length, true));
                int ymid = PartitionY(score_l, score_r);
                Hirschberg(SubCharArray(x,0,xmid,false),SubCharArray(y,0,ymid,false),out zl, out wl);
                Hirschberg(SubCharArray(x,xmid,x.Length,false),SubCharArray(y,ymid,y.Length,false),out zr, out wr);

                z = new char[zl.Length + zr.Length];
                Array.Copy(zl,z,zl.Length);
                Array.Copy(zr,0,z,zl.Length,zr.Length);

                w = new char[wl.Length + wr.Length];
                Array.Copy(wl, w, wl.Length);
                Array.Copy(wr, 0, w, wl.Length, wr.Length);
            }
        }

        /// <summary>
        /// Calculates subarray from char array
        /// </summary>
        /// <param name="x">Input source array</param>
        /// <param name="start">Sampling start index</param>
        /// <param name="end">Sampling end index, exclusive</param>
        /// <param name="shouldReverse">If true, output result is reversed</param>
        /// <returns>The subarray of source char array</returns>
        private static char[] SubCharArray(char[] x, int start, int end, bool shouldReverse)
        {
            char[] result = new char[end-start];
            Array.Copy(x,start,result,0,end-start);
            if(shouldReverse)
                Array.Reverse(result);
            return result;
        }

        /// <summary>
        /// Calculates the alignment score
        /// </summary>
        /// <param name="x">First sequence</param>
        /// <param name="y">Second sequence</param>
        /// <param name="matching">output matching string containing '|' for match, 'x' for mismatch, ' ' for gap</param>
        /// <param name="score">output score</param>
        private static void CalculateScore(char[] x, char[] y, out char[] matching, out int score)
        {
            score = 0;
            matching = new char[x.Length];

            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] == '-' || y[i] == '-')
                {
                    score += GAP;
                    matching[i] = ' ';
                }
                else
                {
                    if (x[i] == y[i])
                    {
                        score += MATCH;
                        matching[i] = '|';
                    }
                    else
                    {
                        score += MISMATCH;
                        matching[i] = 'x';
                    }
                }
            }
        }

        /// <summary>
        /// Alignment process using SmithWaterman, NeedlemanWunsch and Hirschberg algorithm
        /// </summary>
        /// <param name="x">The first sequence</param>
        /// <param name="y">The second sequence</param>
        /// <param name="isWritingToConsole">The flag that, if true, writes steps of alignment to the console, along with time taken</param>
        public static void Align(char[] x, char[] y, bool isWritingToConsole, bool isCheckingMemory)
        {
            char[] xLocal;
            char[] yLocal;
            int maxScore;
            LocalToGlobal(x, y, out xLocal, out yLocal, out maxScore);
            
            if(isCheckingMemory)
                CheckMemory();

            if (isWritingToConsole)
            {
                //Console.WriteLine(maxScore);
                Console.WriteLine("LocalToGlobal:");
                Console.WriteLine(xLocal);
                Console.WriteLine(yLocal);
                Console.WriteLine("-------");
            }
            
            Hirschberg(xLocal, yLocal, out XResult, out YResult);

            if (isCheckingMemory)
                CheckMemory();

            //CalculateScore(XResult, YResult, out Matching, out Score);

            if (isCheckingMemory)
                CheckMemory();

            if (isWritingToConsole)
            {
                //Console.WriteLine("Score: " + Score);
                Console.WriteLine("Alignment result:");
                Console.WriteLine(XResult);
                //Console.WriteLine(Matching);
                Console.WriteLine(YResult);
            }
        }

        /// <summary>
        /// The total recorded amount of memory used by this process
        /// </summary>
        public static long MemoryUsageMax = 0;

        /// <summary>
        /// Checks memory used by this process
        /// </summary>
        public static void CheckMemory()
        {
            long MemoryUsageAcquired = GC.GetTotalMemory(true);
            if (MemoryUsageAcquired > MemoryUsageMax)
                MemoryUsageMax = MemoryUsageAcquired;
        }
     }
}

using System;
using System.IO;
using System.Text;

namespace Hirschberg
{
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
}

using System;

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
    }
}

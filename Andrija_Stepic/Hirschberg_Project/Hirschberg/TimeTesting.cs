using System;
using System.Text;

namespace Hirschberg
{
    public static class TimeTesting
    {
        /// <summary>
        /// Method that tests four ways to copy an array of chararcters for performance
        /// </summary>
        public static void TestCopy()
        {
            Random random = new Random();

            for (int i = 10000; i < 1000000000; i *= 10)
            {
                Console.WriteLine("Started with " + i);
                DateTime startTime = DateTime.Now;

                StringBuilder builder = new StringBuilder();
                for (int j = 0; j < i; j++)
                {
                    builder.Append((char)(random.Next(0, 4) + '0'));
                }

                char[] x = builder.ToString().ToCharArray();

                Console.WriteLine("Time taken to build strings: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");

                char[] y;

                Console.WriteLine("String length: " + x.Length);
                
                startTime = DateTime.Now;
                CopyNaive(x, out y);
                Console.WriteLine("Naive Copy Time taken: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");

                /*
                startTime = DateTime.Now;
                CopyPtr(x, out y);
                Console.WriteLine("Ptr *Copy* Time taken: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");
                */

                startTime = DateTime.Now;
                Copy(x, out y);
                Console.WriteLine("Array.Copy Time taken: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");

                startTime = DateTime.Now;
                CopyBuffer(x, out y,x.Length);
                Console.WriteLine("BufferCopy Time taken: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");
            }
        }

        public static void CopyNaive(char[] src, out char[] dst)
        {
            dst = new char[src.Length];
            for (int i = 0; i < dst.Length; i++)
                dst[i] = src[i];
        }

        public static void Copy(char[] src, out char[] dst)
        {
            dst = new char[src.Length];
            Array.Copy(src, dst, src.Length);
        }
        /*
        public static unsafe void CopyPtr(char[] src, out char[] dst)
        {
            dst = new char[src.Length];
            fixed (char* pSrc = src, pDst = dst)
            {
                char* ps = pSrc;
                char* pd = pDst;
                for (int i = 0; i < src.Length; i++)
                {
                    *pd = *ps;
                    ps++;
                    pd++;
                }
            }
        }
        */
        public static void CopyBuffer(char[] src, out char[] dst, int length)
        {
            dst = new char[src.Length];
            Buffer.BlockCopy(src,0,dst,0,length * sizeof(char));
        }

        public static char[] SubCharArray(char[] x, int start, int end,bool shouldReverse)
        {
            char[] result = new char[end - start];
            Array.Copy(x, start, result, 0, end - start);
            if (shouldReverse)
                Array.Reverse(result);
            return result;
        }

        public static char[] SubCharArrayBuffer(char[] src, int start, int end, bool shouldReverse)
        {
            char[] dst = new char[end - start];
            Buffer.BlockCopy(src, start * sizeof(char), dst, 0, (end - start) * sizeof(char));
            if(shouldReverse)
                Array.Reverse(dst);
            return dst;
        }

        /// <summary>
        /// Method that tests performance of Sequence alignment, using randomly generated
        /// sequences with lengths from 100 to 1000000, both the same length and 10% difference
        /// of characters between them
        /// </summary>
        public static void TestHirschberg()
        {
            Random random = new Random();

            float changePercentage = 0.5f;

            for (int i = 100; i < 1000000; i *= 10)
            {
                Console.WriteLine("Started with " + i);
                DateTime startTime = DateTime.Now;

                StringBuilder builder = new StringBuilder();
                for (int j = 0; j < i; j++)
                {
                    builder.Append((char)(random.Next(0, 4) + '0'));
                }

                char[] x = builder.ToString().ToCharArray();

                for (int j = 0; j < (int)(i * changePercentage); j++)
                {
                    int index = random.Next() % builder.Length;
                    builder[index] = (char)('0' + '3' - builder[index]);
                }

                char[] y = builder.ToString().ToCharArray();

                Console.WriteLine("Time taken to build strings: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");

                SequenceAlignment.Align(x, y, false, false);
                Console.WriteLine("String length: " + Math.Max(x.Length, y.Length));
                Console.WriteLine("Time taken: " + (DateTime.Now - startTime).TotalMilliseconds + "ms");
                //Console.WriteLine("Memory: " + SequenceAlignment.MemoryUsageMax);
            }
        }
    }
}

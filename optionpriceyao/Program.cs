using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OptionPricing2
{
    class Program
    {
        static void Main(string[] args)
        {

            // Get essential class's instances.
            OptionPricing.Output ResultOutput = new OptionPricing.Output();

            //Initial parameters
            double S, K, r, sigma, T, div;
            int N;
            string chose;

            // Create interface
            Console.WriteLine("Trinomial tree calculator");
            Console.WriteLine("");
            Console.WriteLine("Put the value into calculator and pick the option");
            Console.WriteLine("");
            Console.WriteLine("Get the option price and correspond greeks");
            Console.WriteLine("");

            Console.WriteLine("Enter underlying price 'S'");

            S = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Enter Strike 'K'");
            K = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Enter Interest Rate 'r'");
            r = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Enter volatility 'sigma'");
            sigma = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Enter Tenor 'T'");
            T = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Enter Tree Steps 'N'");
            N = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Enter Dividend 'div'");
            div = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine();


            
            Console.WriteLine("Please choose A-D to choose property of the option");
            Console.WriteLine("A.  European Call Option");
            Console.WriteLine("B.  European Put Option");
            Console.WriteLine("C.  American Call Option");
            Console.WriteLine("D.  American Put Option");
            
            chose = Console.ReadLine();


            switch (chose)
            {
                case "A":
                    // Define properties of class Output
                    ResultOutput.IsEuropean = true;
                    ResultOutput.IsCall = true;
                    // Output the result
                    ResultOutput.Output_Calculation(S, K, r, T, sigma, N, div);
                    Console.WriteLine("The calculation is over, thank you for using it.");
                    break;
                case "B":
                    // Define properties of class Output
                    ResultOutput.IsEuropean = true;
                    ResultOutput.IsCall = false;
                    // Output the result
                    ResultOutput.Output_Calculation(S, K, r, T, sigma, N, div);
                    Console.WriteLine("The calculation is over, thank you for using it.");

                    break;
                case "C":
                    // Define properties of class Output
                    ResultOutput.IsEuropean = false;
                    ResultOutput.IsCall = true;
                    // Output the result
                    ResultOutput.Output_Calculation(S, K, r, T, sigma, N, div);
                    Console.WriteLine("The calculation is over, thank you for using it.");

                    break;

                case "D":
                    // Define properties of class Output
                    ResultOutput.IsEuropean = false;
                    ResultOutput.IsCall = false;
                    // Output the result
                    ResultOutput.Output_Calculation(S, K, r, T, sigma, N, div);
                    Console.WriteLine("The calculation is over, thank you for using it.");

                    break;
                // If something wrong happened.
                default:
                    Console.WriteLine("Attention! You didn't follow the instruction, punish you input from beginning!!");
                    break;
            }


        }
    }
}
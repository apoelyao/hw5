using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OptionPricing
{
    /// <summary>
    /// This is a class to get greeks of option, except vega and rho.
    /// </summary>
    class Greek : Option
    {
        // Define greek's properties
        public double[,] Vt { get; set; }
        public double[,] St { get; set; }

        // Get delta
        public double Delta()
        {
            double delta;
            delta = (Vt[0, 1] - Vt[2, 1]) / (St[0, 1] - St[2, 1]);
            delta = Math.Round(delta, 4);
            return delta;
        }

        // Get gamma
        public double Gamma()
        {
            double gamma;
            gamma = (((Vt[0, 1] - Vt[1, 1]) / (St[0, 1] - St[1, 1])) - ((Vt[1, 1] - Vt[2, 1]) / (St[1, 1] - St[2, 1]))) / (0.5 * (St[0, 1] - St[2, 1]));
            gamma = Math.Round(gamma, 4);
            return gamma;
        }

        // Get theta
        public double Theta(double T, int N)
        {
            double theta = (Vt[2, 2] - Vt[0, 0]) / (2 * (T / N));
            theta = Math.Round(theta, 4);
            return theta;
        }
    }

    /// <summary>
    /// This is a class including the method to calculate vega and rho
    /// </summary>
    public class VegaRhoCalculation
    {
        public Boolean Iscall { get; set; }
        public Boolean IsEuropean { get; set; }
        double[,] St_max, St_min, St;

        // Get vega
        public double[,] GetVegaRho(double diff_sigma, double diff_r, double S, double K, double r, double sigma, double T, int N, double div)
        {
            // Initial parameters that used in the following code.
            TrinomialUnderlyingTree tree = new TrinomialUnderlyingTree();
            double vega, rho;
            double[,] vega_rho = new double[1, 2];
            double sigma_max = sigma + diff_sigma;
            double sigma_min = sigma - diff_sigma;
            double r_max = r + diff_r;
            double r_min = r - diff_r;
            // Get underlying price trinomial tree under different sigma.
            St_max = tree.MakeUnderlyingTree(S, sigma_max, T, N);
            St_min = tree.MakeUnderlyingTree(S, sigma_min, T, N);
            St = tree.MakeUnderlyingTree(S, sigma, T, N);

            // Judge if this is an European Option
            if (IsEuropean == true)
            {
                // Define class EurOption instance.
                EurOption eur = new EurOption();
                //Judge if this is a call option.
                if (Iscall == true)
                {
                    eur.Iscall = true;
                    // Calculate Vega
                    vega = (eur.MakingOptionValueTree(K, r, sigma_max, T, N, div, St_max)[0, 0] - eur.MakingOptionValueTree(K, r, sigma_min, T, N, div, St_min)[0, 0]) / (2 * diff_sigma);
                    vega = Math.Round(vega, 4);
                    // Calculate rho
                    rho = (eur.MakingOptionValueTree(K, r_max, sigma, T, N, div, St)[0, 0] - eur.MakingOptionValueTree(K, r_min, sigma, T, N, div, St)[0, 0]) / (2 * diff_r);
                    rho = Math.Round(rho, 4);
                }
                // if this is a put option.
                else
                {
                    eur.Iscall = false;
                    // Calculate Vega
                    vega = (eur.MakingOptionValueTree(K, r, sigma_max, T, N, div, St_max)[0, 0] - eur.MakingOptionValueTree(K, r, sigma_min, T, N, div, St_min)[0, 0]) / (2 * diff_sigma);
                    vega = Math.Round(vega, 4);
                    // Calculate rho
                    rho = (eur.MakingOptionValueTree(K, r_max, sigma, T, N, div, St)[0, 0] - eur.MakingOptionValueTree(K, r_min, sigma, T, N, div, St)[0, 0]) / (2 * diff_r);
                    rho = Math.Round(rho, 4);
                }
            }
            // If this is an American Option
            else
            {
                // Define class AmericanOption instance.
                AmericanOption ame = new AmericanOption();
                //Judge if this is a call option.
                if (Iscall == true)
                {
                    ame.Iscall = true;
                    // Calculate Vega
                    vega = (ame.MakingOptionValueTree(K, r, sigma_max, T, N, div, St_max)[0, 0] - ame.MakingOptionValueTree(K, r, sigma_min, T, N, div, St_min)[0, 0]) / (2 * diff_sigma);
                    vega = Math.Round(vega, 4);
                    // Calculate rho
                    rho = (ame.MakingOptionValueTree(K, r_max, sigma, T, N, div, St)[0, 0] - ame.MakingOptionValueTree(K, r_min, sigma, T, N, div, St)[0, 0]) / (2 * diff_r);
                    rho = Math.Round(rho, 4);
                }
                // if this is a put option.
                else
                {
                    ame.Iscall = false;
                    // Calculate Vega
                    vega = (ame.MakingOptionValueTree(K, r, sigma_max, T, N, div, St_max)[0, 0] - ame.MakingOptionValueTree(K, r, sigma_min, T, N, div, St_min)[0, 0]) / (2 * diff_sigma);
                    vega = Math.Round(vega, 4);
                    // Calculate rho
                    rho = (ame.MakingOptionValueTree(K, r_max, sigma, T, N, div, St)[0, 0] - ame.MakingOptionValueTree(K, r_min, sigma, T, N, div, St)[0, 0]) / (2 * diff_r);
                    rho = Math.Round(rho, 4);
                }
            }
            // Output the result
            vega_rho[0, 0] = vega;
            vega_rho[0, 1] = rho;
            return vega_rho;

        }
    }




    /// <summary>
    /// This is a class to output all calculations of option.
    /// </summary>
    public class Output
    {
        public Boolean IsEuropean { get; set; }
        public Boolean IsCall { get; set; }

        public void Output_Calculation(double S, double K, double r, double T, double sigma, int N, double div)
        {
            // Get essential classes' instances.
            PreCalculate precal = new PreCalculate();
            TrinomialUnderlyingTree tree = new TrinomialUnderlyingTree();
            Option opt = new Option();
            EurOption Eur_option = new EurOption();
            AmericanOption Amer_option = new AmericanOption();
            Greek greeks = new Greek();
            VegaRhoCalculation vega_cal = new VegaRhoCalculation();

            // Initial array St and Vt
            double[,] St = new double[2 * N + 1, N + 1];
            double[,] Vt = new double[2 * N + 1, N + 1];

            // Initial interval that used in Greeks
            double diff_sigma = 0.001 * sigma;
            double diff_r = 0.001 * r;

            // Judge if this is an European Option
            if (IsEuropean == true)
            {
                //Judge if this is a call option.
                if (IsCall == true)
                {
                    Console.WriteLine("The Underlying Prices Trinomial Tree:");
                    // Initial the parameters
                    Eur_option.Iscall = true;
                    vega_cal.IsEuropean = true;
                    vega_cal.Iscall = true;

                    // Ilerate underlying trinomial tree St.
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            St[i, j] = Math.Round(tree.MakeUnderlyingTree(S, sigma, T, N)[i, j], 2);
                            Console.Write(St[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }

                    Console.WriteLine();
                    Console.WriteLine("The Option Value Trinomial Tree:");

                    // Ilerate option value trinomial tree Vt
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            Vt[i, j] = Math.Round(Eur_option.MakingOptionValueTree(K, r, sigma, T, N, div, St)[i, j], 2);
                            Console.Write(Vt[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }
                    // Output the option price 
                    Console.WriteLine("European Call Option Price is" + "  " + Vt[0, 0].ToString());
                    Console.WriteLine();

                    // Output Greeks
                    Console.WriteLine("Greeks of this option ");
                    greeks.St = St;
                    greeks.Vt = Vt;
                    double delta = greeks.Delta();
                    Console.WriteLine("Delta =" + "  " + delta.ToString());
                    double gamma = greeks.Gamma();
                    Console.WriteLine("Gamma =" + "  " + gamma.ToString());
                    double theta = greeks.Theta(T, N);
                    Console.WriteLine("Theta =" + "  " + theta.ToString());
                    double vega = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 0];
                    double rho = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 1];
                    Console.WriteLine("Vega =" + "  " + vega.ToString());
                    Console.WriteLine("Rho =" + "  " + rho.ToString());
                    Console.ReadLine();
                }
                // If this is a put option.
                else
                {
                    Console.WriteLine("The Underlying Prices Trinomial Tree:");
                    // Initial the parameters
                    Eur_option.Iscall = false;
                    vega_cal.IsEuropean = true;
                    vega_cal.Iscall = false;
                    // Ilerate underlying trinomial tree St.
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            St[i, j] = Math.Round(tree.MakeUnderlyingTree(S, sigma, T, N)[i, j], 2);
                            Console.Write(St[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }

                    Console.WriteLine();
                    Console.WriteLine("The Option Value Trinomial Tree:");

                    // Ilerate option value trinomial tree Vt
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            Vt[i, j] = Math.Round(Eur_option.MakingOptionValueTree(K, r, sigma, T, N, div, St)[i, j], 2);
                            Console.Write(Vt[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }
                    // Output option price.
                    Console.WriteLine("European Put Option Price is" + "  " + Vt[0, 0].ToString());
                    Console.WriteLine();

                    // Output Greeks
                    Console.WriteLine("Greeks of this option :");
                    greeks.St = St;
                    greeks.Vt = Vt;
                    double delta = greeks.Delta();
                    Console.WriteLine("Delta =" + "  " + delta.ToString());
                    double gamma = greeks.Gamma();
                    Console.WriteLine("Gamma =" + "  " + gamma.ToString());
                    double theta = greeks.Theta(T, N);
                    Console.WriteLine("Theta =" + "  " + theta.ToString());
                    double vega = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 0];
                    double rho = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 1];
                    Console.WriteLine("Vega =" + "  " + vega.ToString());
                    Console.WriteLine("Rho =" + "  " + rho.ToString());
                    Console.ReadLine();
                }
            }
            // If this is an American Option
            else
            {
                //Judge if this is a call option.
                if (IsCall == true)
                {
                    Console.WriteLine("The Underlying Prices Trinomial Tree:");
                    // Initial the parameters
                    Amer_option.Iscall = true;
                    vega_cal.IsEuropean = false;
                    vega_cal.Iscall = true;
                    // Ilerate underlying trinomial tree St.
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            St[i, j] = Math.Round(tree.MakeUnderlyingTree(S, sigma, T, N)[i, j], 2);
                            Console.Write(St[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }

                    Console.WriteLine();
                    Console.WriteLine("The Option Value Trinomial Tree:");

                    // Ilerate option value trinomial tree Vt
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            Vt[i, j] = Math.Round(Amer_option.MakingOptionValueTree(K, r, sigma, T, N, div, St)[i, j], 2);
                            Console.Write(Vt[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }
                    // Output American call option price
                    Console.WriteLine("American Call Option Price is" + "  " + Vt[0, 0].ToString());
                    Console.WriteLine();

                    // Output Greeks
                    Console.WriteLine("Greeks of this option:");
                    greeks.St = St;
                    greeks.Vt = Vt;
                    double delta = greeks.Delta();
                    Console.WriteLine("Delta =" + "  " + delta.ToString());
                    double gamma = greeks.Gamma();
                    Console.WriteLine("Gamma =" + "  " + gamma.ToString());
                    double theta = greeks.Theta(T, N);
                    Console.WriteLine("Theta =" + "  " + theta.ToString());
                    double vega = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 0];
                    double rho = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 1];
                    Console.WriteLine("Vega =" + "  " + vega.ToString());
                    Console.WriteLine("Rho =" + "  " + rho.ToString());
                    Console.ReadLine();
                }
                //If this is a put option.
                else
                {
                    Console.WriteLine("The Underlying Prices Trinomial Tree:");
                    // Initial the parameters
                    Amer_option.Iscall = false;
                    vega_cal.IsEuropean = false;
                    vega_cal.Iscall = false;
                    // Ilerate underlying trinomial tree St.
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            St[i, j] = Math.Round(tree.MakeUnderlyingTree(S, sigma, T, N)[i, j], 2);
                            Console.Write(St[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }

                    Console.WriteLine();
                    Console.WriteLine("The Option Value Trinomial Tree:");

                    // Ilerate option value trinomial tree Vt
                    for (int i = 0; i <= 2 * N; i++)
                    {

                        for (int j = 0; j <= N; j++)
                        {
                            Vt[i, j] = Math.Round(Amer_option.MakingOptionValueTree(K, r, sigma, T, N, div, St)[i, j], 2);
                            Console.Write(Vt[i, j].ToString() + "\t");

                        }
                        Console.WriteLine();
                    }

                    // Output American call option price
                    Console.WriteLine("American Put Option Price is" + "  " + Vt[0, 0].ToString());
                    Console.WriteLine();

                    // Output Greeks
                    Console.WriteLine("Greeks of this option:");
                    greeks.St = St;
                    greeks.Vt = Vt;
                    double delta = greeks.Delta();
                    Console.WriteLine("Delta =" + "  " + delta.ToString());
                    double gamma = greeks.Gamma();
                    Console.WriteLine("Gamma =" + "  " + gamma.ToString());
                    double theta = greeks.Theta(T, N);
                    Console.WriteLine("Theta =" + "  " + theta.ToString());
                    double vega = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 0];
                    double rho = vega_cal.GetVegaRho(diff_sigma, diff_r, S, K, r, sigma, T, N, div)[0, 1];
                    Console.WriteLine("Vega =" + "  " + vega.ToString());
                    Console.WriteLine("Rho =" + "  " + rho.ToString());
                    Console.ReadLine();
                }
            }
        }
    }
}



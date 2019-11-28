using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OptionPricing
{
    /// <summary>
    /// This is a class to calculate serval parameters that will use in the future.
    /// </summary>
    class PreCalculate
    {
        public double sigma { get; set; }
        public double r { get; set; }
        public double T { get; set; }
        public int N { get; set; }
        public double div { get; set; }

        // Get time interval dt
        public double Get_dt()
        {
            return (T / N);
        }

        // Get drift fatctor nu
        public double Get_nu()
        {
            return (r - div - 0.5 * Math.Pow(sigma, 2));
        }

        // Get size of moving price dx
        public double Get_dx()
        {
            double dt;
            dt = Get_dt();
            return (sigma * Math.Sqrt(3 * dt));
        }

        // Get continuous moving price edx
        public double Get_edx()
        {
            double dx;
            dx = Get_dx();
            return (Math.Exp(dx));
        }

        // Get discount rate
        public double Get_disc()
        {
            double dt;
            dt = Get_dt();
            return Math.Exp(-r * dt);

        }

        // pu, pm, pd are probabilities under risk-neutral measure.
        // Pre-calculate the pu, pm, pd of discount to improve efficiency of the loop in the following code.
        //Get pu after discount
        public double Get_disc_pu()
        {
            double pu, dt, dx, nu, disc, disc_pu;
            dt = Get_dt();
            nu = Get_nu();
            dx = Get_dx();
            disc = Get_disc();

            pu = 0.5 * (((Math.Pow(sigma, 2) * dt + Math.Pow(nu, 2) * Math.Pow(dt, 2)) / Math.Pow(dx, 2)) + ((nu * dt) / dx));
            disc_pu = disc * pu;
            return disc_pu;
        }

        // Get pm after discount
        public double Get_disc_pm()
        {
            double pm, dt, dx, nu, disc, disc_pm;
            dt = Get_dt();
            nu = Get_nu();
            dx = Get_dx();
            disc = Get_disc();

            pm = 1.0 - ((Math.Pow(sigma, 2) * dt + Math.Pow(nu, 2) * Math.Pow(dt, 2)) / Math.Pow(dx, 2));
            disc_pm = disc * pm;
            return disc_pm;
        }

        // Get pd after discount
        public double Get_disc_pd()
        {
            double pd, dt, dx, nu, disc, disc_pd;
            dt = Get_dt();
            nu = Get_nu();
            dx = Get_dx();
            disc = Get_disc();

            pd = 0.5 * ((Math.Pow(sigma, 2) * dt + Math.Pow(nu, 2) * Math.Pow(dt, 2)) / Math.Pow(dx, 2) - (nu * dt) / dx);
            disc_pd = disc * pd;
            return disc_pd;
        }


    }

    /// <summary>
    /// This is a class to get underlying price trinomial tree.
    /// </summary>
    public class TrinomialUnderlyingTree
    {
        public double[,] MakeUnderlyingTree(double S, double sigma, double T, int N)
        {
            PreCalculate precal = new PreCalculate();
            precal.sigma = sigma;
            precal.T = T;
            precal.N = N;

            // Get St on every node.
            // Get discount rate
            double edx = precal.Get_edx();
            // Get size of moving price dx
            double dx = precal.Get_dx();
            //Initiate the outcome St.
            double[,] St = new double[2 * N + 1, N + 1];

            // Calculate the underlying price if it always goes down in the whole period.
            St[2 * N, N] = S * Math.Exp(-N * dx);

            // Calculate all the possible underlying price at the last step.
            for (int j = (2 * N - 1); j >= 0; j--)
            {
                St[j, N] = St[j + 1, N] * edx;
            }

            // Calculate all underlying prices St.
            for (int i = N - 1; i >= 0; i--)
            {
                for (int j = 2 * i; j >= 0; j--)
                {
                    // Every node, except that on the last step, has node on the next step which has the same price.
                    St[j, i] = St[j + 1, i + 1];
                }
            }

            return St;
        }


    }

    /// <summary>
    /// This is a class to define option properties
    /// </summary>
    class Option
    {
        public Boolean Iscall { get; set; }

        public virtual double[,] MakingOptionValueTree(double K, double r, double sigma, double T, int N, double div, double[,] StTree)
        {
            double[,] OptionValue = new double[1, 1];
            return OptionValue;
        }
    }

    /// <summary>
    /// This is a class to get European call and put option prices.
    /// </summary>
    class EurOption : Option
    {
        public override double[,] MakingOptionValueTree(double K, double r, double sigma, double T, int N, double div, double[,] StTree)
        {
            // Initial class Precalculate
            PreCalculate precal = new PreCalculate();
            precal.sigma = sigma;
            precal.r = r;
            precal.T = T;
            precal.N = N;
            precal.div = div;
            // Get pu, pm, pd under discount
            double disc_pu = precal.Get_disc_pu();
            double disc_pm = precal.Get_disc_pm();
            double disc_pd = precal.Get_disc_pd();

            //Initiate St and outcome Vt.
            double[,] St = new double[2 * N + 1, N + 1];
            double[,] Vt = new double[2 * N + 1, N + 1];
            // Copy old array that got from method to the new array.
            for (int i = 0; i <= 2 * N; i++)
            {

                for (int j = 0; j <= N; j++)
                {
                    St[i, j] = StTree[i, j];


                }
            }

            // Calculate the option value of last step if it is a call option.
            if (Iscall == true)
            {
                for (int j = (2 * N); j >= 0; j--)
                {
                    Vt[j, N] = Math.Max(St[j, N] - K, 0);
                }

            }
            // Calculate the option value of last step if it is a put option.
            else
            {
                for (int j = (2 * N); j >= 0; j--)
                {
                    Vt[j, N] = Math.Max(K - St[j, N], 0);
                }
            }


            for (int i = N - 1; i >= 0; i--)
            {
                for (int j = 2 * i; j >= 0; j--)
                {
                    // Every node, except that on the last step, equals expectation of nodes' on last step.
                    Vt[j, i] = disc_pu * Vt[j, i + 1] + disc_pm * Vt[j + 1, i + 1] + disc_pd * Vt[j + 2, i + 1];
                }
            }

            return Vt;

        }

    }

    /// <summary>
    /// This is a class to get American call and put option prices.
    /// </summary>
    class AmericanOption : Option
    {
        public override double[,] MakingOptionValueTree(double K, double r, double sigma, double T, int N, double div, double[,] StTree)
        {
            // Initial class Precalculate
            PreCalculate precal = new PreCalculate();
            precal.sigma = sigma;
            precal.r = r;
            precal.T = T;
            precal.N = N;
            precal.div = div;
            // Get pu, pm, pd under discount
            double disc_pu = precal.Get_disc_pu();
            double disc_pm = precal.Get_disc_pm();
            double disc_pd = precal.Get_disc_pd();

            //Initiate St and outcome Vt.
            double[,] St = new double[2 * N + 1, N + 1];
            double[,] Vt = new double[2 * N + 1, N + 1];
            // Copy old array that got from method to the new array.
            for (int i = 0; i <= 2 * N; i++)
            {

                for (int j = 0; j <= N; j++)
                {
                    St[i, j] = StTree[i, j];


                }
            }

            // Calculate the option value of last step if it is a call option.
            if (Iscall == true)
            {
                for (int j = (2 * N); j >= 0; j--)
                {
                    Vt[j, N] = Math.Max(St[j, N] - K, 0);
                }

                // Calculate the option value of every nodes.
                for (int i = N - 1; i >= 0; i--)
                {
                    for (int j = 2 * i; j >= 0; j--)
                    {
                        // Every node, except that on the last step, equals expectation of nodes' on last step.
                        Vt[j, i] = disc_pu * Vt[j, i + 1] + disc_pm * Vt[j + 1, i + 1] + disc_pd * Vt[j + 2, i + 1];
                        // Compare the value.
                        Vt[j, i] = Math.Max(Vt[j, i], St[j, i] - K);
                    }
                }
            }

            // Calculate the option value of last step if it is a put option.
            else
            {
                for (int j = (2 * N); j >= 0; j--)
                {
                    Vt[j, N] = Math.Max(K - St[j, N], 0);
                }

                // Calculate the option value of every nodes.
                for (int i = N - 1; i >= 0; i--)
                {
                    for (int j = 2 * i; j >= 0; j--)
                    {
                        // Every node, except that on the last step, equals expectation of nodes' on last step.
                        Vt[j, i] = disc_pu * Vt[j, i + 1] + disc_pm * Vt[j + 1, i + 1] + disc_pd * Vt[j + 2, i + 1];
                        // Compare the value.
                        Vt[j, i] = Math.Max(Vt[j, i], K - St[j, i]);
                    }
                }
            }
            return Vt;

        }

    }

}
package sceua;

public class Problem {

    public double problem(int nopt, double[] x, int problemIndex){

        Double result = null;

        if(problemIndex == 1){
            result = goldsteinPrice(nopt, x);
        }else if (problemIndex == 2) {
            result = rosenbrock(nopt, x);
        }else if (problemIndex == 3) {
            result = sixHumpCamelback(nopt, x);
        }else if (problemIndex == 4) {
            result = rastrigin(nopt, x);
        }else if (problemIndex == 5) {
            result = griewank(nopt, x);
        }else if (problemIndex == 6) {
            result = shekel(nopt, x);
        }else if (problemIndex == 7) {
            result = hartman(nopt, x);
        }else {
            System.out.println("未找到代码对应的问题");
        }

        return result;
    }

    // This is the Goldstein-Price Function
    // Bound X1=[-2,2], X2=[-2,2]
    // Global Optimum: 3.0,(0.0,-1.0)
    public double goldsteinPrice(int nopt, double[] x){

        double x1 = x[0];
        double x2 = x[1];
        double u1 = Math.pow(x1 + x2 + 1.0, 2);
        double u2 = 19.0 - 14 * x1 + 3.0 * Math.pow(x1, 2)
                - 14.0 * x2 + 6.0 * x1 * x2 + 3.0 * Math.pow(x2, 2);
        double u3 = Math.pow(2.0*x1 - 3.0*x2, 2);
        double u4 = 18.0 - 32.0*x1 + 12.0*Math.pow(x1, 2) + 48.0*x2
                - 36.0*x1*x2 + 27.0*Math.pow(x2, 2);
        double u5 = u1 * u2;
        double u6 = u3 * u4;
        double f = (1.0 + u5) * (30.0 + u6);

        return f;
    }

    //%  This is the Rosenbrock Function
    //%  Bound: X1=[-5,5], X2=[-2,8]
    //%  Global Optimum: 0,(1,1)
    public double rosenbrock(int nopt, double[] x){

        double x1 = x[0];
        double x2 = x[1];
        double a = 100;
        double b = 1;
        double f = a * Math.pow(x2 - Math.pow(x1, 2), 2) + Math.pow((b - x1), 2);

        return f;
    }

    //%  This is the Six-hump Camelback Function.
    //%  Bound: X1=[-5,5], X2=[-5,5]
    //%  True Optima: -1.031628453489877, (-0.08983,0.7126), (0.08983,-0.7126)
    public double sixHumpCamelback(int nopt, double[] x){

        double x1 = x[0];
        double x2 = x[1];
        double term1 = 4 - 2.1 * Math.pow(x1, 2) + Math.pow(x1, 4) / 3;
        double term2 = x1 * x1;
        double term3 = x1 * x2;
        double term4 = -4 + 4 * Math.pow(x2, 2);
        double term5 = x2 * x2;

        double f = term1 * term2 + term3 + term4 * term5;

        return f;

    }

    //%  This is the Rastrigin Function
    //%  Bound: X1=[-1,1], X2=[-1,1]
    //%  Global Optimum: -2, (0,0)
    public double rastrigin(int nopt, double[] x){

        double x1 = x[0];
        double x2 = x[1];
        double term1 = Math.pow(x1, 2);
        double term2 = Math.pow(x2, 2);
        double term3 = Math.cos(18.0 * x1);
        double term4 = Math.cos(18.0 * x2);

        double f = term1 + term2 - term3 - term4;

        return f;
    }

    //%  This is the Griewank Function (2-D or 10-D)
    //%  Bound: X(i)=[-600,600], for i=1,2,...,10
    //%  Global Optimum: 0, at origin
    public double griewank(int nopt, double[] x){

        int d;

        if(nopt == 2){
            d = 200;
        }else {
            d = 4000;
        }

        double u1 = 0.0;
        double u2 = 1.0;

        for (int i = 0; i < nopt; i++) {
            u1 = u1 + Math.pow(x[i], 2) / d;
            u2 = u2 * Math.cos(x[i] / Math.sqrt(i + 1)); // Note: MATLAB uses 1-based indexing
        }

        double f = u1 - u2 + 1;

        return f;
    }

    //%  This is the Shekel Function
    //%  Bound: X(i)=[0,10], j=1,2,3,4
    //%  Global Optimum:-10.5364098252,(4,4,4,4)
    public double shekel(int nopt, double[] x){

        // Data for Skekel function coefficients (n=4, m=10)
        double[][] a1 = {
                {4., 1., 8., 6., 3., 2., 5., 8., 6., 7.},
                {4., 1., 8., 6., 7., 9., 5., 1., 2., 3.6},
                {4., 1., 8., 6., 3., 2., 3., 8., 6., 7.},
                {4., 1., 8., 6., 7., 9., 3., 1., 2., 3.6}
        };

        double[] c1 = {.1, .2, .2, .4, .4, .6, .3, .7, .5, .5};

        double f = 0.0;

        for (int i = 0; i < 10; i++) {
            double u = 0.0;

            for (int j = 0; j < nopt; j++) {
                u = u + Math.pow(x[j] - a1[j][i], 2);
            }

            u = 1.0 / (u + c1[i]);
            f = f - u;
        }

        return f;
    }

    //%   This is the Hartman Function
    //%   Bound: X(j)=[0,1], j=1,2,...,6
    //%   Global Optimum:-3.322368011415515,
    //%   (0.201,0.150,0.477,0.275,0.311,0.657)
    public double hartman(int nopt, double[] x){

        // Data for Hartman function coefficients (6-D)
        double[][] a2 = {
                {10.0, 0.05, 3.0, 17.0},
                {3.0, 10.0, 3.5, 8.0},
                {17.0, 17.0, 1.7, 0.05},
                {3.5, 0.1, 10.0, 10.0},
                {1.7, 8.0, 17.0, 0.1},
                {8.0, 14.0, 8.0, 14.0}
        };
        double[] c2 = {1.0, 1.2, 3.0, 3.2};
        double[][] p2 = {
                {0.1312, 0.2329, 0.2348, 0.4047},
                {0.1696, 0.4135, 0.1451, 0.8828},
                {0.5569, 0.8307, 0.3522, 0.8732},
                {0.0124, 0.3736, 0.2883, 0.5743},
                {0.8283, 0.1004, 0.3047, 0.1091},
                {0.5886, 0.9991, 0.6650, 0.0381}
        };
        // Data for Hartman function coefficient (3-D)
        double[][] a3 = {
                {3.0, 0.1, 3.0, 0.1},
                {10.0, 10.0, 10.0, 10.0},
                {30.0, 35.0, 30.0, 35.0}
        };
        double[] c3 = {1.0, 1.2, 3.0, 3.2};
        double[][] p3 = {
                {0.3689, 0.4699, 0.1091, 0.03815},
                {0.1170, 0.4387, 0.8732, 0.5743},
                {0.2673, 0.7470, 0.5547, 0.8828}
        };

        double f = 0.0;

        for (int i = 0; i < 4; i++) {
            double u = 0.0;

            for (int j = 0; j < nopt; j++) {
                if (nopt == 3) {
                    a2[j][i] = a3[j][i];
                    p2[j][i] = p3[j][i];
                }
                u = u + a2[j][i] * Math.pow(x[j] - p2[j][i], 2);
            }

            if (nopt == 3) {
                c2[i] = c3[i];
            }

            f = f - c2[i] * Math.exp(-u);
        }

        return f;

    }



}

package sceua;

import org.apache.commons.math3.random.RandomDataGenerator;

public class CCEUA {

    /**  This is the subroutine for generating a new point in a simplex
    //%   s(.,.) = the sorted simplex in order of increasing function values
    //%   s(.) = function values in increasing order
    //%  LIST OF LOCAL VARIABLES
    //%   sb(.) = the best point of the simplex
    //%   sw(.) = the worst point of the simplex
    //%   w2(.) = the second worst point of the simplex
    //%   fw = function value of the worst point
    //%   ce(.) = the centroid of the simplex excluding worst
    //%   snew(.) = new point generated from the simplex
    //%   iviol = flag indicating if constraints are violated
    //%         = 1 , yes
    //%         = 0 , no
     */
    public void cceua(double[][] s, double[] sf, double[] bl, double[] bu,
                      int maxn,// 没用上的值
                      double[] snew, double[] fnew, int[] icall,// 三个返回值
                      int problemIndex,// 目标函数代码
                      RandomDataGenerator randomDataGenerator){

        int nps = s.length;
        int nopt = s[0].length;
        int n = nps;
        int m = nopt;
        double alpha = 1.0;
        double beta = 0.5;

        //% Assign the best and worst points:
        double[] sb = s[0];double fb = sf[0];
        double[] sw = s[n-1];double fw = sf[n-1];

        //% Compute the centroid of the simplex excluding the worst point:
        // ce=mean(s(1:n-1,:));
        MatlabFunction matf = new MatlabFunction();
        // 先复制前n-1行
        double[][] sTmp = new double[n-1][];
        for (int i = 0; i < n - 1; i++) {
            sTmp[i] = s[i];
        }
        double[] ce = matf.matlabMean(sTmp);

        //% Attempt a reflection point
//        double[] snew = new double[ce.length];
        // snew在外面初始化 直接操作snew
        for (int i = 0; i < m; i++) {
            snew[i] = ce[i] + alpha * (ce[i] - sw[i]);
        }

        //% Check if is outside the bounds:
        int ibound = 0;
        double[] s1 = new double[snew.length];
        for (int i = 0; i < snew.length; i++) {
            s1[i] = snew[i] - bl[i];
        }
        if(hasNegativeElements(s1)){
            ibound=1;
        }
        // 更新S1内的值
        for (int i = 0; i < snew.length; i++) {
            s1[i] = bu[i] - snew[i];
        }
        if(hasNegativeElements(s1)){
            ibound=2;
        }

        // 有负数则重新生成一个解
        if (ibound >= 1) {
            for (int i = 0; i < snew.length; i++) {
                snew[i] = randomDataGenerator.nextUniform(bl[i], bu[i]);
            }
        }

        Problem pb = new Problem();

        // 操作fnew
        fnew[0] = pb.problem(nopt, snew, problemIndex);
        icall[0] = icall[0] + 1;

        //% Reflection failed; now attempt a contraction point:
        if (fnew[0] > fw) {// 比最大loss还大
            for (int i = 0; i < m; i++) {
                snew[i] = sw[i] + beta * (ce[i] - sw[i]);
            }
            fnew[0] = pb.problem(nopt, snew, problemIndex);
            icall[0] = icall[0] + 1;

            //% Both reflection and contraction have failed, attempt a random point;
            if (fnew[0] > fw) {// 还比最差的差
                for (int i = 0; i < ce.length; i++) {
                    snew[i] = randomDataGenerator.nextUniform(bl[i], bu[i]);
                }
                fnew[0] = pb.problem(nopt, snew, problemIndex);
                icall[0] = icall[0]+ 1;
            }
        }

    }

    // 检查数组是否存在小于零的元素
    private static boolean hasNegativeElements(double[] arr) {
        for (double value : arr) {
            if (value < 0) {
                return true; // 如果找到小于零的元素，返回 true
            }
        }
        return false; // 如果数组中所有元素都大于或等于零，返回 false
    }
}

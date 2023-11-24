import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.ArrayList;
import java.util.Arrays;


public void sceua(double[] x0, double[] bl, double[] bu,
                      int maxn, int kstop, double pcento, double peps,
                      int ngs, int iseed, int iniflg,
                      int problemIndex){// 目标函数代码

        // Initialize SCE parameters:
        int nopt = x0.length;
        int npg = 2 * nopt + 1;
        int nps = nopt + 1;
        int nspl = npg;
//        int mings = ngs;
        int npt = npg * ngs;// 所有解的数量

        double[] bound = new double[nopt];
        for (int i = 0; i < bound.length; i++) {
            bound[i] = bu[i] - bl[i];
        }

        // Create an initial population to fill array x(npt,nopt):
        RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        randomDataGenerator.reSeed(1);  // 设置种子值
        double[][] x = new double[npt][nopt];
        for (int i = 0; i < npt; i++) {
            for (int j = 0; j < x[i].length; j++) {
                // 生成随机数默认不包含上下界
                x[i][j] = randomDataGenerator.nextUniform(bl[j], bu[j]);
            }
        }

        // 是否包含初始解
        if (iniflg == 1){
            x[0] = x0;
        }

        Problem pb = new Problem();

        int nloop = 0;
        int icall = 0;
        double[] xf = new double[npt];
        for (int i = 0; i < npt; i++) {
//            xf[i] = function();// 每个解对应的目标函数值
            xf[i] = pb.problem(nopt, x[i], problemIndex);
            icall = icall + 1;
        }
//        f0=xf(1);



        // Sort the population in order of increasing function values;
        Integer[] idx = new Integer[x.length];
        // 初始化 idx 数组
        for (int i = 0; i < x.length; i++) {
            idx[i] = i;
        }
        // 使用 xf 对 idx 进行排序 idx是排序前元素位置
        matf.matlabSort(xf, idx);
        // 根据idx对x排序
        double[][] sortedX = new double[x.length][];
        for (int i = 0; i < x.length; i++) {
            sortedX[i] = x[idx[i]];
        }
        // 替换原先的x
        x = sortedX;

        // Record the best and worst points;
        double[] bestx = x[0];// xf最小值
        double bestf = xf[0];
        double[] worstx = x[npt-1];
        double worstf = xf[npt-1];
        BESTF.add(bestf);
        Double[] wrapperBestx = new Double[bestx.length];
        for (int i = 0; i < bestx.length; i++) {
            wrapperBestx[i] = bestx[i]; // 自动装箱
        }
        BESTX.add(wrapperBestx);
        ICALL.add(icall);

        // Compute the standard deviation for each parameter
//        StandardDeviation stdDev = new StandardDeviation();
//        double[] xnstd = new double[nopt];
//        for (int i = 0; i < xnstd.length; i++) {
//            xnstd[i] = stdDev.evaluate(x[i]);
//        }
        double[] xnstd = matf.matlabStd(x);


        // Computes the normalized geometric range of the parameters
//        gnrng=exp(mean(log((max(x)-min(x))./bound)));
        // 计算每一列的最小值和最大值
        double[] minValues = new double[nopt];// 解空间的维度 解向量的长度
        double[] maxValues = new double[nopt];
        for (int col = 0; col < nopt; col++) {
            // 初始化为第一行的值
            minValues[col] = x[0][col];
            maxValues[col] = x[0][col];
            // 遍历每一行
            for (int row = 1; row < x.length; row++) {
                double value = x[row][col];
                // 更新最小值和最大值
                minValues[col] = Math.min(minValues[col], value);
                maxValues[col] = Math.max(maxValues[col], value);
            }
        }
        // 计算 (max(x) - min(x)) ./ bound
        double[] ratios = new double[nopt];
        for (int i = 0; i < ratios.length; i++) {
            ratios[i] = (maxValues[i] - minValues[i]) / bound[i];
        }
        // 计算 log((max(x) - min(x)) / bound)
        double[] logResults = new double[ratios.length];
        for (int i = 0; i < ratios.length; i++) {
            logResults[i] = Math.log(ratios[i]);
        }
        // 计算 exp(mean(log((max(x) - min(x)) / bound)))
        double mean = matf.matlabMean(logResults);
        double gnrng = Math.exp(mean);

        System.out.println("The Initial Loop: 0");
        System.out.println("BESTF  : " + bestf);
        System.out.println("BESTX  : " + Arrays.toString(bestx));
        System.out.println("WORSTF : " + worstf);
        System.out.println("WORSTX : " + Arrays.toString(worstx));
        System.out.println("------------------------------");

        // Check for convergency;
        if (icall >= maxn) {
            System.out.println("*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT");
            System.out.println("ON THE MAXIMUM NUMBER OF TRIALS ");
            System.out.println(maxn);
            System.out.println("HAS BEEN EXCEEDED.  SEARCH WAS STOPPED AT TRIAL NUMBER:");
            System.out.println(icall);
            System.out.println("OF THE INITIAL LOOP!");
        }

        if (gnrng < peps){
            System.out.println("THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE");
        }

        // Begin evolution loops:
        nloop = 0;
        ArrayList<Double> criter = new ArrayList<>();
        double criter_change = 1e+5;

        while (icall<maxn & gnrng>peps & criter_change>pcento){

            nloop = nloop + 1;

            // Loop on complexes (sub-populations);
            for(int igs = 1; igs <= ngs; igs++){

                // Partition the population into complexes (sub-populations);
                int[] k1 = new int[npg];
                for (int i = 0; i < npg; i++) {
                    k1[i] = i + 1;
                }
                int[] k2 = new int[k1.length];
                for (int i = 0; i < k1.length; i++) {
                    k2[i] = (k1[i] - 1) * ngs + igs;
                }
                // 创建并初始化 cx 和 cf 数组
                double[][] cx = new double[npg][nopt];
                double[] cf = new double[npg];
                // 复制数组元素
                for (int i = 0; i < k1.length; i++) {// K1K2长度相同
                    for (int j = 0; j < nopt; j++) {
                        cx[k1[i]-1][j] = x[k2[i]-1][j];
                    }
                    cf[k1[i]-1] = xf[k2[i]-1];
                }
                // Evolve sub-population igs for nspl steps:
                for(int loop = 1; loop <= nspl; loop++){// nspl=npg
                    // Select simplex by sampling the complex according to a linear
                    // probability distribution
                    int[] lcs = new int[nps];
                    lcs[0] = 1;
                    for(int k3 = 2; k3 <= nps; k3++) {// k3和i一样
                        int lpos = 1;
                        for (int iter = 1; iter <= 1000; iter++) {
                            lpos = 1 +  (int) Math.floor(
                                    npg + 0.5 - Math.sqrt(
                                            Math.pow(npg + 0.5, 2) - npg * (npg + 1) * randomDataGenerator.nextUniform(0, 1)
                                    )
                            );
                            // 保证之前的元素中没有本次生成的lpos
                            int idx2 = find(lcs, lpos, k3);
                            if (idx2 == -1) {
                                break;
                            }
                        }
                        lcs[k3 - 1] = lpos;
                    }
                    Arrays.sort(lcs);

                    //% Construct the simplex:
                    //            s=cx(lcs,:); sf = cf(lcs);
                    double[][] s = new double[nps][nopt];
                    double[] sf = new double[nps];
                    for (int i = 0; i < nps; i++) {
                        // 将 cx 数组中由 lcs[i] 索引确定的行复制到 s 数组中
                        for (int j = 0; j < nopt; j++) {// cx长npg s和lcs长nps
                            s[i][j] = cx[lcs[i]-1][j];
                        }
                        // 将 cf 数组中由 lcs[i] 索引确定的元素复制到 sf 数组中
                        sf[i] = cf[lcs[i] - 1];
                    }

                    // 调用CCEUA
                    double[] snew = new double[nopt];
                    double[] fnew = new double[1];// 包装一下用方法传值
                    int[] wrappedIcall = {icall};// 包装一下用方法传值
                    cceua.cceua(s, sf, bl, bu, maxn,
                            snew, fnew, wrappedIcall, // snew, fnew, icall为返回值
                            problemIndex, randomDataGenerator);
                    icall = wrappedIcall[0];

                    //% Replace the worst point in Simplex with the new point:
                    s[nps-1] = snew;
                    sf[nps-1] = fnew[0];

                    //% Replace the simplex into the complex;
                    //  cx(lcs,:) = s;
                    //  cf(lcs) = sf;
                    // 将 s和sf 的值赋给 cx和cf 的指定行
                    for (int i = 0; i < lcs.length; i++) {// cx长npg s和lcs长nps
                        int rowIndex = lcs[i] - 1;
                        // 检查 rowIndex 是否在合法范围内
                        if (rowIndex >= 0 && rowIndex < cx.length) {
                            // 更新 cx 的指定行
                            cx[rowIndex] = Arrays.copyOf(s[i], s[i].length);
                            cf[rowIndex] = sf[i];
                        } else {
                            System.out.println("Invalid row index: " + rowIndex);
                        }
                    }

                    //% Sort the complex;
                    // [cf,idx] = sort(cf); cx=cx(idx,:);
                    Integer[] idx3 = new Integer[cf.length];
                    // 初始化 idx3 数组
                    for (int i = 0; i < cf.length; i++) {
                        idx3[i] = i;
                    }
                    // 使用 cf 对 idx3 进行排序
                    matf.matlabSort(cf, idx3);
                    // 根据idx对x排序
                    double[][] sortedCx = new double[cf.length][];
                    for (int i = 0; i < cf.length; i++) {
                        sortedCx[i] = cx[idx3[i]];
                    }
                    // 替换原先的x
                    cx = sortedCx;
                    //% End of Inner Loop for Competitive Evolution of Simplexes
                }

                //% Replace the complex back into the population;
                // x(k2,:) = cx(k1,:);
                // xf(k2) = cf(k1);
                // 复制数组元素 k1k2长度相同
                for (int i = 0; i < k1.length; i++) {
                    for (int j = 0; j < nopt; j++) {
                        x[k2[i]-1][j] = cx[k1[i]-1][j];
                    }
                    xf[k2[i]-1] = cf[k1[i]-1];
                }

                //% End of Loop on Complex Evolution;
            }

            //% Shuffled the complexes;
            //    [xf,idx] = sort(xf); x=x(idx,:);
            //    PX=x; PF=xf;
            idx = new Integer[xf.length];
            // 初始化 idx 数组
            for (int i = 0; i < x.length; i++) {
                idx[i] = i;
            }
            // 使用 xf 对 idx 进行排序 idx是排序前元素位置
            matf.matlabSort(xf, idx);
            // 根据idx对x排序
            sortedX = new double[x.length][];
            for (int i = 0; i < x.length; i++) {
                sortedX[i] = x[idx[i]];
            }
            // 替换原先的x
            x = sortedX;

            PX=x;
            PF=xf;

            //% Record the best and worst points;
            //    bestx=x(1,:); bestf=xf(1);
            //    worstx=x(npt,:); worstf=xf(npt);
            bestx = x[0];// xf最小值
            bestf = xf[0];
            worstx = x[npt-1];
            worstf = xf[npt-1];
            //    BESTX=[BESTX;bestx]; BESTF=[BESTF;bestf];ICALL=[ICALL;icall];
            BESTF.add(bestf);
            wrapperBestx = new Double[bestx.length];// 重新装填
            for (int i = 0; i < bestx.length; i++) {
                wrapperBestx[i] = bestx[i]; // 自动装箱
            }
            BESTX.add(wrapperBestx);
            ICALL.add(icall);

            //% Compute the standard deviation for each parameter
            //xnstd=std(x);
            xnstd = matf.matlabStd(x);

            //% Computes the normalized geometric range of the parameters
            //    gnrng=exp(mean(log((max(x)-min(x))./bound)));
            minValues = new double[nopt];// 解空间的维度 解向量的长度
            maxValues = new double[nopt];
            for (int col = 0; col < nopt; col++) {
                // 初始化为第一行的值
                minValues[col] = x[0][col];
                maxValues[col] = x[0][col];
                // 遍历每一行
                for (int row = 1; row < x.length; row++) {
                    double value = x[row][col];
                    // 更新最小值和最大值
                    minValues[col] = Math.min(minValues[col], value);
                    maxValues[col] = Math.max(maxValues[col], value);
                }
            }
            // 计算 (max(x) - min(x)) ./ bound
            ratios = new double[nopt];
            for (int i = 0; i < ratios.length; i++) {
                ratios[i] = (maxValues[i] - minValues[i]) / bound[i];
            }
            // 计算 log((max(x) - min(x)) / bound)
            logResults = new double[ratios.length];
            for (int i = 0; i < ratios.length; i++) {
                logResults[i] = Math.log(ratios[i]);
            }
            // 计算 exp(mean(log((max(x) - min(x)) / bound)))
            mean = matf.matlabMean(logResults);
            gnrng = Math.exp(mean);

            System.out.println("Evolution Loop: " + nloop + "  - Trial - " + icall);
            System.out.println("BESTF  : " + bestf);
            System.out.println("BESTX  : " + Arrays.toString(bestx));
            System.out.println("WORSTF : " + worstf);
            System.out.println("WORSTX : " + Arrays.toString(worstx));
            System.out.println("------------------------------");

            //% Check for convergency;
            //    if icall >= maxn;
            //        disp('*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT');
            //        disp(['ON THE MAXIMUM NUMBER OF TRIALS ' num2str(maxn) ' HAS BEEN EXCEEDED!']);
            //    end;
            if (icall >= maxn) {
                System.out.println("*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE LIMIT");
                System.out.println("ON THE MAXIMUM NUMBER OF TRIALS " + maxn + " HAS BEEN EXCEEDED!");
            }

            if (gnrng < peps) {
                System.out.println("THE POPULATION HAS CONVERGED TO A PRESPECIFIED SMALL PARAMETER SPACE");
            }

            criter.add(bestf);
            //if (nloop >= kstop);
            //        criter_change=abs(criter(nloop)-criter(nloop-kstop+1))*100;
            //        criter_change=criter_change/mean(abs(criter(nloop-kstop+1:nloop)));
            //        if criter_change < pcento;
            //            disp(['THE BEST POINT HAS IMPROVED IN LAST ' num2str(kstop) ' LOOPS BY ', ...
            //                  'LESS THAN THE THRESHOLD ' num2str(pcento) '%']);
            //            disp('CONVERGENCY HAS ACHIEVED BASED ON OBJECTIVE FUNCTION CRITERIA!!!')
            //        end;
            //    end;
            if (nloop >= kstop){
                criter_change = Math.abs(criter.get(nloop-1) - criter.get(nloop - kstop)) * 100;
                //abs(criter(nloop-kstop+1:nloop)
                double sum = 0.0;
                for (int i = nloop - kstop ; i < nloop; i++) {
                    sum = sum + Math.abs(criter.get(i));
                }
                //mean(abs(criter(nloop-kstop+1:nloop)))
                double criMean = sum / kstop;
                //criter_change/mean(abs(criter(nloop-kstop+1:nloop)));
                criter_change = criter_change / criMean;

                if (criter_change < pcento) {
                    System.out.println("THE BEST POINT HAS IMPROVED IN LAST " + kstop + " LOOPS BY "
                            + "LESS THAN THE THRESHOLD " + pcento + "%");
                    System.out.println("CONVERGENCY HAS ACHIEVED BASED ON OBJECTIVE FUNCTION CRITERIA!!!");
                }
            }

            //% End of the Outer Loops
        }

        System.out.println("SEARCH WAS STOPPED AT TRIAL NUMBER: " + icall);
        System.out.println("NORMALIZED GEOMETRIC RANGE = " + gnrng);
        System.out.println("THE BEST POINT HAS IMPROVED IN LAST " + kstop + " LOOPS BY "
                + criter_change + "%");
        //% END of Subroutine sceua
    }


    // 在数组的前 k3-1 个元素中查找是否存在指定的值
    private static int find(int[] arr, int value, int endIdx) {
        for (int i = 0; i < endIdx; i++) {
            if (arr[i] == value) {
                return i;
            }
        }
        return -1;
    }

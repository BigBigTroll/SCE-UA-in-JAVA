package sceua;

import java.util.Arrays;
import java.util.Comparator;

public class MatlabFunction {
    public double[] matlabStd(double[][] dataIn){

        double[] columnStdDevs = new double[dataIn[0].length];

        for (int col = 0; col < dataIn[0].length; col++) {
            double[] colValue = new double[dataIn.length];
            for (int row = 0; row < dataIn.length; row++) {
                colValue[row] = dataIn[row][col];
            }
            // mean
            double sum = 0.0;
            for (double value : colValue) {
                sum = sum + value;
            }
            double mean = sum / colValue.length;
            // 计算差值的平方和
            double sumSquaredDiff = 0.0;
            for (double value : colValue) {
                double diff = value - mean;
                sumSquaredDiff = sumSquaredDiff + diff * diff;
            }
            // 计算方差
            double variance = sumSquaredDiff / (colValue.length - 1);
            // 计算标准差
            columnStdDevs[col] = Math.sqrt(variance);
        }

        return columnStdDevs;
    }

    public double[] matlabMean(double[][] dataIn){

        double[] result = new double[dataIn[0].length];
        // 遍历每一列
        for (int col = 0; col < dataIn[0].length; col++) {
            double sum = 0;
            // 遍历前 n-1 行
            for (int row = 0; row < dataIn.length; row++) {
                sum += dataIn[row][col];
            }
            // 计算平均值，存储在 ce 数组中
            result[col] = sum / dataIn.length;
        }

        return result;
    }

    public double matlabMean(double[] dataIn){

        double sum = 0;

        for (int i = 0; i < dataIn.length; i++) {
            sum = sum + dataIn[i];
        }

        return sum / dataIn.length;
    }

    public void matlabSort(double[] dataIn, Integer[] idx){

        // 使用 dataIn 对 idx 进行排序
        Arrays.sort(idx, new Comparator<Integer>() {
            @Override
            public int compare(final Integer o1, final Integer o2) {
                return Double.compare(dataIn[o1], dataIn[o2]);
            }
        });

        // 新数组排序
        Arrays.sort(dataIn);
    }


  

}

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Main {
    /**
     * 特征矩阵
     */
    private double[][] feature;
    /**
     * 样本标签
     */
    private int[] label;
    /**
     * 梯度下降法步长
     */
    private double stepLength;
    /**
     * 最大迭代次数
     */
    private int maxStep;
    /**
     * 权重矩阵初始化值
     */
    private double initWeight;
    /**
     * 训练后的权重矩阵
     */
    private double[] weights;

    public double[] getWeights() {
        return weights;
    }
    // 训练数据
    private String trainFileName;
    // 测试数据
    private String testFileName;
    // 预测结果
    private String predictFileName;

    public Main(String trainFileName, String testFileName, String predictFileName) {
        this.trainFileName = trainFileName;
        this.testFileName = testFileName;
        this.predictFileName = predictFileName;

        this.stepLength = 0.1;
        this.maxStep = 300;
        this.initWeight = 1.0;
    }

    private void loadTrainingData() {
        double[][] matrix = loadFile(trainFileName, false);

        feature = new double[matrix.length][matrix[0].length];
        label = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length - 1; j++) {
                feature[i][j] = matrix[i][j];
            }
            label[i] = (int) matrix[i][matrix[i].length - 1];
        }
    }

    private void initWeightMatrix() {
        int paraSize = feature[0].length;
        double[] weights = new double[paraSize];
        for (int i = 0; i < paraSize; i++) {
            weights[i] = initWeight;
        }
        this.weights = weights;
    }

    /**
     * 迭代标签值
     *
     * @return 预测标签值
     */
    private double[] getPredictLabel() {
        double[] predictLabels = new double[feature.length];
        for (int i = 0; i < predictLabels.length; i++) {
            double predictSum = 0;
            for (int j = 0; j < feature[i].length; j++) {
                predictSum += feature[i][j] * weights[j];
            }
            predictLabels[i] = sigmoid(predictSum);
        }
        return predictLabels;
    }

    /**
     * 计算权重矩阵偏差
     *
     * @return 权重矩阵偏差
     */
    private double[] getDeltaWeights() {
        double[] predictLabels = getPredictLabel();
        double[] deltaWeights = new double[feature[0].length];
        for (int i = 0; i < feature[0].length; i++) {
            deltaWeights[i] = 0;
            for (int j = 0; j < feature.length; j++) {
                deltaWeights[i] += feature[j][i] * (label[j] - predictLabels[j]);
            }
            deltaWeights[i] /= feature.length;
        }

        // System.out.println("Error: "+getErrorRate(predictLabels));

        return deltaWeights;
    }

    // 计算误差率
    private double getErrorRate(double[] predictLabels) {
        double sumErr = 0.0;
        for (int i = 0; i < label.length; i++) {
            sumErr += Math.pow(this.label[i] - predictLabels[i], 2);
        }
        return sumErr;
    }

    /**
     * 迭代训练权重矩阵
     */
    public void training() {
        loadTrainingData();

        if (feature.length <= 0 || feature[0].length <= 0) {
            weights = null;
            return;
        }
        initWeightMatrix();
        for (int i = 0; i < maxStep; i++) {
            double[] deltaWeight = getDeltaWeights();
            for (int j = 0; j < feature[0].length; j++) {
                weights[j] += stepLength * deltaWeight[j];
            }
        }
        this.feature = null;
        this.label = null;
    }

    public double[][] loadFile(String fileName, boolean skipTitle) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(fileName));
        } catch (FileNotFoundException exception) {
            System.err.println(fileName + " File Not Found");
            return null;
        }
        List<List<Double>> listArr = new ArrayList<>();
        String line = "";
        try {
            if (skipTitle) {
                reader.readLine();
            }
            while ((line = reader.readLine()) != null) {
                List<Double> list = new ArrayList<>();
                String item[] = line.split(",");
                for (int i = 0; i < item.length; i++) {
                    list.add(Double.parseDouble(item[i]));
                }
                listArr.add(list);
            }
        } catch (IOException exception) {
            System.err.println(exception.getMessage());
        }

        double[][] matrix = new double[listArr.size()][listArr.get(0).size()];
        for (int i = 0; i < listArr.size(); i++) {
            for (int j = 0; j < listArr.get(i).size(); j++) {
                matrix[i][j] = listArr.get(i).get(j);
            }
        }
        return matrix;
    }

    public void predict() {
        double[][] testFeature = loadFile(testFileName, false);
        int[] predictLabel = new int[testFeature.length];
        for (int i = 0; i < testFeature.length; i++) {
            double sum = 0;
            for (int j = 0; j < testFeature[i].length; j++) {
                sum += testFeature[i][j] * weights[j];
            }
            predictLabel[i] = sigmoid(sum) > 0.5 ? 1 : 0;
        }
        savePredictResult(predictLabel);
    }

    private void savePredictResult(int[] predictLabel) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(predictFileName));
            for (int i = 0; i < predictLabel.length; i++) {
                out.write(predictLabel[i] + "\n");
            }
            out.close();
        } catch (IOException exception) {
            System.err.println(exception.getMessage());
        }
    }

    private double sigmoid(double x) {
        return 1.0d / (1.0d + Math.exp(-x));
    }

    // 命令行格式为test train_data.txt（训练集） test_data.txt(测试数据) predict.txt(预测结果) [debug]
    public static void main(String[] args) {
        BufferedReader reader = null;

        String trainFileName = "/data/train_data.txt";
        String testFileName = "/data/test_data.txt";
        String predictFileName = "/projects/student/result.txt";
        String answerFileName = "/projects/student/answer.txt";

        boolean isDebug = false;
        if (args.length >= 1) {
            if (args[0].equals("debug")) {
                isDebug = true;
            }
        }

        long start = System.currentTimeMillis();
        Main lr = new Main(trainFileName, testFileName, predictFileName);
        long end = System.currentTimeMillis();
        System.out.println("Init Time(s): " + (end - start) * 1.0 / 1000);
        lr.training();
        System.out.println("Training Time(s): " + (System.currentTimeMillis() - end) * 1.0 / 1000);
        end = System.currentTimeMillis();
        lr.predict();
        System.out.println("Predict Time(s): " + (System.currentTimeMillis() - end) * 1.0 / 1000);
        end = System.currentTimeMillis();
        if (isDebug) {
            double[][] matrix = lr.loadFile(predictFileName, false);
            int[] predict = new int[matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                predict[i] = (int) matrix[i][0];
            }
            matrix = lr.loadFile(answerFileName, false);
            int[] answer = new int[matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                answer[i] = (int) matrix[i][0];
            }

            int accCount = 0;
            for (int i = 0; i < predict.length; i++) {
                if (predict[i] == answer[i]) {
                    accCount++;
                }
            }
            System.out.println("Accuracy:" + (1.0f * accCount / predict.length));
            System.out.println("Mark Time(s): " + (System.currentTimeMillis() - end) * 1.0 / 1000);
        }
    }
}

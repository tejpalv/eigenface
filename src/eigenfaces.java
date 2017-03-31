import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class eigenfaces {
	// image name 1.pgm to 7.pgm
	static int start = 1;
	static int end = 7;
	// dimensions of image:
	static int width = 140; 
	static int height = 198;
	
	public static void main(String[] args) {
		
		double[][] faces = new double[end][width * height]; //matrix of faces in the form of column vectors of n^2 x 1
		for(int i = start; i <= end; i++){
				PgmImage img;
				String filename ="/Users/tejpalvirdi/Desktop/eigenfae/" + i + ".pgm";
				img = new PgmImage(filename);
				double[] columnVector = toColumnVector(img.returnPixels());
				for(int j = 0; j < columnVector.length; j++){
					faces[i-1][j] = columnVector[j];
				}
	}
		double[] average = average(faces); //computes average 
				
		double[][] twoDdoubles = columnTo2D(average (faces));
		
		
		double[][] facesMinusAverage = facesMinusAverage(faces, average);		
		
		double[][] innerProduct = getInnerProduct(facesMinusAverage);
		
			Matrix INNERPRODUCT = new Matrix(innerProduct);
			
		    EigenvalueDecomposition eig = INNERPRODUCT.eig();
		    Matrix A = new Matrix(facesMinusAverage); 
		    A = A.transpose(); // how i implemented it puts A in the wrong orientation (transpose)
			double[][] eigenvectors = eig.getV().getArray();
			double[][] oneVector = new double[eigenvectors.length][eigenvectors.length];
			
				for(int i = 0; i < eigenvectors.length; i++){
					for(int j=0; j<eigenvectors.length; j++){
						oneVector[j][i] = eigenvectors[j][i]; // moving eigenvectors to single row vector
					}
				}

			Matrix singleVECTOR = new Matrix(eigenvectors);
			singleVECTOR = singleVECTOR.transpose(); // oneVector transposed
		
			Matrix tempVECTOR = A.times(singleVECTOR); // eigenvectors multiplied by tranpose

			double[][] fff = tempVECTOR.getArray();
			int[][] zzz = new int[fff.length][fff[0].length];
			for(int i = 0; i < fff.length; i++){
				for(int j = 0; j <fff[0].length; j++)
				{
					zzz[i][j] = (int) ((int) 255 * fff[i][j]); //eigenfactors scaled to 255 and converted to ints;
				}
			}
			
			tempVECTOR = tempVECTOR.transpose();
			
			double[][] weights = new double[faces.length][faces.length];
						
			double[][] prom = new double[1][facesMinusAverage[0].length];
			
			
			// Ω = U^T(Γ−Ψ)
			for(int i=0 ; i<faces.length ; i++){
				prom[0] = facesMinusAverage[3];
				Matrix pro = new Matrix(prom);
				pro = pro.transpose();
				weights = tempVECTOR.times(pro).getArray();
				
			}

		
			 //reconstruct:
			for(int i=0; i<faces.length;i++){
				
				for(int j=0;j<faces[0].length;j++){
					faces[i][j] = weights[i][0] * faces[i][j];
				}
			}
			int[][] p = convertToIntArr(faces);
			
			// empty image:
			int[][] blank = new int[27720][1]; 
						
			// RENDERING:
			
			//reset:
			//writePGM(blank);
				
			//Weights:
			System.out.println(Arrays.deepToString(weights));

			//Average Face:
			//writePGM(convertToIntArr(twoDdoubles));
			
			//Prominent Parts Face:
			//writePGM(convertToIntArr(facesMinusAverage));
			
			//Eigenface:
			//writePGM(zzz);
			
			//Reconstruction:
			//writePGM(p);
			
	}

private static double[][] getInnerProduct(double[][] facesMinusAverage) {
		// TODO Auto-generated method stub

	double[][] innerProduct = new double[facesMinusAverage.length][facesMinusAverage.length];
	for(int i = 0; i < innerProduct.length; i++){
		for(int j = 0; j < innerProduct[0].length; j++){
			int sum = 0;
			for(int k = 0; k < width * height; k++){
				sum += facesMinusAverage[i][k] * facesMinusAverage[j][k];
			}
			innerProduct[i][j] = sum;
		}
	}
	return innerProduct;
	}


	private static int[][] convertToIntArr(double[][] twoDdoubles) {
		// TODO Auto-generated method stub
		int[][] intarr = new int[twoDdoubles.length][twoDdoubles[0].length];
		
		for(int i=0; i<twoDdoubles.length;i++){
			for(int j=0; j<twoDdoubles[0].length;j++){
				intarr[i][j] = (int) twoDdoubles[i][j];
			}
		}
		
		return intarr;
	}


	public static double[] toColumnVector(int[][] x){
		double[] z = new double[x.length * x[0].length];
		for(int j = 0; j < x.length; j++){
			for(int i = 0; i < x[0].length; i++){
				z[(j * x[0].length) + i] = x[j][i];
			}
		}
		return z;
		
	}
	
	public static double[] average(double[][] faces){
		double[] w = new double[faces[0].length];
		for(int i = 0; i < faces[0].length; i++){
			int tempAVG = 0;
			for(int j = 0; j < faces.length; j++){
				tempAVG += faces[j][i];
			}
			tempAVG = tempAVG / faces.length;
			w[i] = tempAVG-30;
		}
		return w;
	}
	
	public static double[][] columnTo2D(double[] vector){
		double[][] zeta = new double[198][140];
		for(int j = 0; j <198; j++){
			for(int i = 0; i < 140; i++){
				zeta[j][i] = vector[j * 140 + i];
			}
		}
		return zeta;
	}
	
	public static double[][] facesMinusAverage(double[][] faces, double[] average){
		double[][] facesMinusAverage = faces;
		for(int i = start; i <= end; i++){
			for(int j = 0; j < (width * height); j++){
				facesMinusAverage[i-1][j] =  facesMinusAverage[i-1][j] - (average[j]) -20;
						;
				if(facesMinusAverage[i-1][j] < 0){
					facesMinusAverage[i-1][j] = 0;
				}
			}
		}
		return facesMinusAverage;
	}
	
	public static void writePGM(int[][] ds){
		try{
		     FileWriter fstream = new FileWriter("test.pgm");

		     BufferedWriter out = new BufferedWriter(fstream);
		
		     out.write("P2 \n140 198\n255\n");
		     
		     for(int i = 0 ; i< ds.length;i++){
		        for(int j = 0 ; j<ds[0].length;j++){
		            out.write(ds[i][j] + " ");
		            out.flush();
		            }
		            
		        }
		     
		     }
		catch (Exception e){
		     System.err.println("Error : " + e.getMessage());
		}
		
	}
	
	}

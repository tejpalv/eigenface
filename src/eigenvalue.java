import java.util.Arrays;


public class eigenvalue {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int[][] identity = createIdentity(2,1);
		
		int[][] matrix = {{1,3},{4,5}};
		
		int j=1;
		while(determinant(subtract(matrix , createIdentity(2,j))) != 0){
			j++;
		}
		int eigenvalue = j;
		System.out.println(eigenvalue);
	}

	public static int determinant(int[][] matrix){ //method sig. takes a matrix (two dimensional array), returns determinant.
	    int sum=0; 
	    int s;
	    if(matrix.length==1){  //bottom case of recursion. size 1 matrix determinant is itself.
	      return(matrix[0][0]);
	    }
	    for(int i=0;i<matrix.length;i++){ //finds determinant using row-by-row expansion
	      int[][]smaller= new int[matrix.length-1][matrix.length-1]; //creates smaller matrix- values not in same row, column
	      for(int a=1;a<matrix.length;a++){
	        for(int b=0;b<matrix.length;b++){
	          if(b<i){
	            smaller[a-1][b]=matrix[a][b];
	          }
	          else if(b>i){
	            smaller[a-1][b-1]=matrix[a][b];
	          }
	        }
	      }
	      if(i%2==0){ //sign changes based on i
	        s=1;
	      }
	      else{
	        s=-1;
	      }
	      sum+=s*matrix[0][i]*(determinant(smaller)); 
	      //recursive step: determinant of larger determined by smaller.
	    }
	    return(sum); //returns determinant value. once stack is finished, returns final determinant.
	  }

	public static int[][] subtract(int[][] a, int[][] b) {
	       int rows = a.length;
	       int columns = a[0].length;
	       int[][] result = new int[rows][columns];
	       for (int i = 0; i < rows; i++) {
	           for (int j = 0; j < columns; j++) {
	               result[i][j] = a[i][j] - b[i][j];
	           }
	       }
	       return result;
	   }
	
	private static int[][] createIdentity(int a, int b) {
		int[][] matrix = new int[a][a];
		for(int i=0; i<a;i++){
			for(int j=0; j<a ; j++){
				if(i==j){
					matrix[i][j] = b;
				} else{
					matrix[i][j] = 0;
				}
			}
		}
		
		
		return matrix;
		// TODO Auto-generated method stub
		
	}

}

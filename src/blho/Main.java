package blho;

import java.io.File;
import java.util.Scanner;

public class Main {
	
	
	public static int w[];
	public static double H[][], P[][], PH[][];
	public static double M[][];
	public static double wH[][];
	public static double A[][];
	public static double b[][];
	public static double c[][];
	
	
	public static int m, n = 3;

	
	public static double[][] multiplyByMatrix(double[][] m1, double[][] m2) {
        int m1ColLength = m1[0].length; // m1 columns length
        int m2RowLength = m2.length;    // m2 rows length
        if(m1ColLength != m2RowLength) return null; // matrix multiplication is not possible
        int mRRowLength = m1.length;    // m result rows length
        int mRColLength = m2[0].length; // m result columns length
        double[][] mResult = new double[mRRowLength][mRColLength];
        for(int i = 0; i < mRRowLength; i++) {         // rows from m1
            for(int j = 0; j < mRColLength; j++) {     // columns from m2
                for(int k = 0; k < m1ColLength; k++) { // columns from m1
                    mResult[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return mResult;
    }


	public static double[][] invert(double a[][]) 
	{
		int n = a.length;
		double x[][] = new double[n][n];
		double b[][] = new double[n][n];
		int index[] = new int[n];
		for (int i=0; i<n; ++i) 
			b[i][i] = 1;

		// Transform the matrix into an upper triangle
		gaussian(a, index);

		// Update the matrix b[i][j] with the ratios stored
		for (int i=0; i<n-1; ++i)
			for (int j=i+1; j<n; ++j)
				for (int k=0; k<n; ++k)
					b[index[j]][k]
							-= a[index[j]][i]*b[index[i]][k];

		// Perform backward substitutions
		for (int i=0; i<n; ++i) 
		{
			x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
			for (int j=n-2; j>=0; --j) 
			{
				x[j][i] = b[index[j]][i];
				for (int k=j+1; k<n; ++k) 
				{
					x[j][i] -= a[index[j]][k]*x[k][i];
				}
				x[j][i] /= a[index[j]][j];
			}
		}
		return x;
	}

	// Method to carry out the partial-pivoting Gaussian
	// elimination.  Here index[] stores pivoting order.

	public static void gaussian(double a[][], int index[]) 
	{
		int n = index.length;
		double c[] = new double[n];

		// Initialize the index
		for (int i=0; i<n; ++i) 
			index[i] = i;

		// Find the rescaling factors, one from each row
		for (int i=0; i<n; ++i) 
		{
			double c1 = 0;
			for (int j=0; j<n; ++j) 
			{
				double c0 = Math.abs(a[i][j]);
				if (c0 > c1) c1 = c0;
			}
			c[i] = c1;
		}

		// Search the pivoting element from each column
		int k = 0;
		for (int j=0; j<n-1; ++j) 
		{
			double pi1 = 0;
			for (int i=j; i<n; ++i) 
			{
				double pi0 = Math.abs(a[index[i]][j]);
				pi0 /= c[index[i]];
				if (pi0 > pi1) 
				{
					pi1 = pi0;
					k = i;
				}
			}

			// Interchange rows according to the pivoting order
			int itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;
			for (int i=j+1; i<n; ++i) 	
			{
				double pj = a[index[i]][j]/a[index[j]][j];

				// Record pivoting ratios below the diagonal
				a[index[i]][j] = pj;

				// Modify other elements accordingly
				for (int l=j+1; l<n; ++l)
					a[index[i]][l] -= pj*a[index[j]][l];
			}
		}
	}
	public static double getDeterminant(double[][] matrix){ //method sig. takes a matrix (two dimensional array), returns determinant.
		int sum=0; 
		int s;
		if(matrix.length==1){  //bottom case of recursion. size 1 matrix determinant is itself.
			return(matrix[0][0]);
		}
		for(int i=0;i<matrix.length;i++){ //finds determinant using row-by-row expansion
			double[][]smaller = new double[matrix.length-1][matrix.length-1]; //creates smaller matrix- values not in same row, column
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
			sum+=s*matrix[0][i]*(getDeterminant(smaller)); //recursive step: determinant of larger determined by smaller.
		}
		return(sum); //returns determinant value. once stack is finished, returns final determinant.
	}
	
	public static double getDeterminant2(double A[][]){
	        double res = A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[1][0]*A[2][1]*A[0][2]
	        - A[0][2]*A[1][1]*A[2][0] - A[0][1]*A[1][0]*A[2][2] - A[1][2]*A[2][1]*A[0][0];
	        return res;
	}
	 
	public static void show(double matrix[][]){
		for (int i=0;i<matrix.length ;i++) {
			for (int j=0;j<matrix[i].length;j++) {
				System.out.printf("%3d ", (int)matrix[i][j]);
			}
			System.out.println();
		}
	       
	}
	
	public static void show(int matrix[]){
		for (int i=0;i<matrix.length ;i++) {
			System.out.printf("%s ", matrix[i]);
		}
		System.out.println();
	       
	}
	public static void main (String argv[]) throws Exception {
        Scanner in = new Scanner(new File("input.txt"));

        m = 6;
        n = 3;
        int L = 3;
        
        w = new int[m];
        M = new double[L][L];
        H = new double[L][L];
        P = new double[L][L];
        PH = new double[L][L];
        wH = new double[L][L];
        A = new double[L][L];
        b = new double[1][n];
        c = new double[1][n];
        
        double tmp[][] = null;
        
        for (int i=0;i<m;i++) {
                w[i] = in.nextInt();
        }
 
        System.out.println("Input matrix w:");
        show(w);
        
        
        for (int i=0;i<n;i++) {
                for (int j=0;j<n;j++) {
                        H[i][j] = w[(i+j)%m];
                }
        }
       
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++) {
                    wH[i][j] = w[(i+j + 1)%m];
            }
        }
 
        for (int i=0;i<n;i++) {
                for (int j=0;j<n;j++) {
                        if ( i == j ) {
                        	P[i][j] = 1.0;
                        } else if ( i >= j ) {
                                P[i][j] = ((int)(Math.random() * 10)) % 2 + 1;
                                if ( Math.random() > 0.5f ) P[i][j] *= -1.0;
                        }
                }
        }
        
        P = new double[][] {{1, 0, 0}, {-1, 1, 0}, {-1, 0, 1}};
        
        PH = multiplyByMatrix(P, H);
 
//        for (int i=0;i<n;i++){
//                for (int j=0;j<n;j++){
//                        for (int k=0;k<n;k++){
//                                PH[i][j] += P[i][k] * H[k][j];
//                        }
//                }
//        }
 
        double det = getDeterminant(PH);
 
        System.out.println("Det PH:" + det);
        
        M = invert(PH);
        
//        for (int i=0;i<n;i++){
//                for (int j=0;j<n;j++){
//                        M[i][j] = PH[i][j]/det;                
//                }
//        }
        
        double t[][] = multiplyByMatrix(P, wH);
        A = multiplyByMatrix(t, M);
        
//        for (int i=0;i<n;i++){
//            for (int j=0;j<n;j++){
//                    for (int k=0;k<n;k++){
//                            t[i][j] += P[i][k] * wH[k][j];
//                    }
//            }
//        }

//        for (int i=0;i<n;i++){
//            for (int j=0;j<n;j++){
//                    for (int k=0;k<n;k++){
//                            A[i][j] += t[i][k] * M[k][j];
//                    }
//            }
//        }
        
        b = new double[][]{{M[0][0]}, {M[1][0]}, {M[2][0]}};
        c = new double[][]{{P[0][0], P[0][1], P[0][2]}};
        
        System.out.println("H maxtrix:");
        show(H);
        
        System.out.println("P maxtrix:");
        show(P);
        
        System.out.println("M maxtrix:");
        show(M);
 
        
        tmp = multiplyByMatrix(P, H);
        tmp = multiplyByMatrix(tmp, M);

        
        System.out.println("I maxtrix:");
        show(tmp);

        System.out.println("wH maxtrix:");
        show(wH);

        System.out.println("A maxtrix:");
        show(A);

        System.out.println("b maxtrix:");
        show(b);

        System.out.println("c maxtrix:");
        show(c);

        int nth = 7;
        
        double omega[] = new double[nth];
        
        double Atmp[][] = A;
        
        
        tmp = multiplyByMatrix(c, b); 
        
        omega[0] = tmp[0][0];
        
        for (int i=1;i<nth;i++) {
        	if (i > 1) {
        		Atmp = multiplyByMatrix(Atmp, A);
        	}
        	
        	tmp = multiplyByMatrix(c, Atmp);
        	tmp = multiplyByMatrix(tmp, b);
        	
        	omega[i] = tmp[0][0];
        }
        
        
        for (int i=0;i<omega.length;i++) {
        	System.out.printf("omega[%d]=%d\n", i+1, (int)omega[i]);
        }
	}
	
}

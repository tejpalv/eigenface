import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.*;
import java.util.*;
/**
 * This class handles simple processing for PGM images
 * The existing codes can load and display a gray-scale image in PGM format
 * It also contains an implemented image processing method which 
 *   flips the image horizontally.
 * For homework #1, you should only add new code inside this class
 *   without touching any of the existing code.
 * The new code includes: 
 *   - Four new methods for 4 other simple image processing routines.
 *   - New code in the main method to conduct user specified image manipulations
 *  
 *  @author jxue
 *  @version 0.2011-1-24
 *  
 *  Please follow the homework submission guideline to finish the assignment.
 */
public class PgmImage extends Component {
	// image buffer for graphical display
	private BufferedImage img;
	// image buffer for plain gray-scale pixel values
	private int[][] pixels;
	
	// translating raw gray scale pixel values to buffered image for display
	private void pix2img(){
		int g;
		img = new BufferedImage( pixels[0].length, pixels.length, BufferedImage.TYPE_INT_ARGB );
		// copy the pixels values
		for(int row=0; row<pixels.length; ++row)
			for(int col=0; col<pixels[row].length; ++col){
				g = pixels[row][col];
				img.setRGB(col, row, ((255<<24) | (g << 16) | (g <<8) | g));		
			}
	}

	// default constructor with a 3 by 4 image
	public PgmImage(){
		int[][] defaultPixels = {{0, 1, 2, 3}, {4, 5, 6, 7},{8, 9, 10, 11}};
		pixels = new int[defaultPixels.length][defaultPixels[0].length];
		for(int row=0; row < pixels.length; ++row)
			for(int col=0; col < pixels[0].length; ++col)
				pixels[row][col] = (int)(defaultPixels[row][col] * 255.0/12);
		pix2img();
	}
	
	// constructor that loads pgm image from a file
	public PgmImage(String filename) {
		pixels = null;
		readPGM(filename);
		if (pixels != null)
			pix2img();
	}

	// load gray scale pixel values from a PGM format image
	public void readPGM(String filename){
		try {                        		    
		    Scanner infile = new Scanner(new FileReader(filename));
		    // process the top 4 header lines
		    String filetype=infile.nextLine();
		    if (!filetype.equalsIgnoreCase("p2")) {
		    	System.out.println("[readPGM]Cannot load the image type of "+filetype);
		    	return;
		    }
	   	   	infile.nextLine();	   	   	   
	   	   	int cols = 140;
	   	   	int rows = 198;
	   	   	int maxValue = infile.nextInt();	        
	   	   	pixels = new int[rows][cols];	   	       
	   	   	System.out.println("Reading in image from " + filename + " of size " + rows + " by " + cols);
	   	   	// process the rest lines that hold the actual pixel values
	   	   	for (int r=0; r<rows; r++) 
	   	   		for (int c=0; c<cols; c++)
	   	   			pixels[r][c] = (int)(infile.nextInt()*255.0/maxValue);
	   	   	infile.close();
	    } catch(FileNotFoundException fe) {
	    	System.out.println("Had a problem opening a file.");
	    } catch (Exception e) {
	    	System.out.println(e.toString() + " caught in readPPM.");
	    	e.printStackTrace();
	    }
	}
	// overrides the paint method of Component class
	public void paint(Graphics g) {
		// simply draw the buffered image
		g.drawImage(img, 0, 0, getWidth(), getHeight(), this);
	}
	// overrides the method in Component class, to determine the window size
	public Dimension getPreferredSize() {
		if (img == null) {
			return new Dimension(100, 100);
		} else {			
			// make sure the window is not two small to be seen
			return new Dimension(Math.max(100, img.getWidth(null)), 
						Math.max(100, img.getHeight(null)));
		}
	}
	
	/**
	 * flipH: Flip the image in horizontal direction
	 */
	public void flipH(){
		int tmp, sym;
		for(int row=0; row < pixels.length; ++row){
			for(int col=0; col < pixels[row].length/2; ++col){
				// find the column index of the horizontally symmetric pixel
				sym = pixels[row].length-1 - col;
				// swap the pixel value between the two
				tmp = pixels[row][col];
				pixels[row][col] = pixels[row][sym];
				pixels[row][sym] = tmp;
			}
		}
		this.pix2img();
	}

	/**
	 * flipV: Flip the image in vertical direction
	 */
	/***************************************************/
	 // Add your code here
	/***************************************************/

	/**
	 * inverse: Inverse the image pixels (black->white and vice versa)
	 */
	/***************************************************/
	 // Add your code here
	/***************************************************/

	/**
	 * gridV: add white vertical bars to the image
	 * @param lineWidth the width of the vertical bar lines
	 * @param lineSpacing the spacing between adjacent vertical bar lines
	 */
	public void gridV(int lineWidth, int lineSpacing){
		/***************************************************/
		 // Add your code here
		/***************************************************/		
	}
	
	public int[][] returnPixels(){
		return pixels;
	}
}
	
	/**
	 * gridH: add horizontal bars to the image
	 * @param lineWidth the width of the horizontal bar lines
	 * @param lineSpacing the spacing between adjacent horizontal bar lines
	 * @param lineGrayscale the gray scale value of the horizontal bar lines
	 */
	/***************************************************/
	 // Add your code here

	/***************************************************/

	// The main method that will load and process a pgm image, and display the result.

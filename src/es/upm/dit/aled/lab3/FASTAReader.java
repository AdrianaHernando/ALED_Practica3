package es.upm.dit.aled.lab3;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Reads a FASTA file containing genetic information and allows for the search
 * of specific patterns within these data. The information is stored as an array
 * of bytes that contain nucleotides in the FASTA format. Since this array is
 * usually created before knowing how many characters in the origin FASTA file
 * are valid, an int indicating how many bytes of the array are valid is also
 * stored. All valid characters will be at the beginning of the array.
 * 
 * @author mmiguel, rgarciacarmona
 *
 */
public class FASTAReader {

	protected byte[] content;
	protected int validBytes;

	/**
	 * Creates a new FASTAReader from a FASTA file.
	 * 
	 * @param fileName The name of the FASTA file.
	 */
	public FASTAReader(String fileName) {
		try {
			this.readFile(fileName);
		} catch (IOException e) {
			System.out.println(e.getMessage());
			return;
		}
	}

	/*
	 * Helper method to read from a file. It populates the data array with upper
	 * case version of all the nucleotids found in the file. Throws an IOException
	 * if there is a problem accessing the file or the file is to big to fit in an
	 * array.
	 */
	private void readFile(String fileName) throws IOException {
		File f = new File(fileName);
		FileInputStream fis = new FileInputStream(f);
		DataInput fid = new DataInputStream(fis);
		long len = (int) fis.getChannel().size();
		if (len > Integer.MAX_VALUE) {
			fis.close();
			throw new IOException("The file " + fileName + " is too big. Can't be contained in an array.");
		}
		byte[] content = new byte[(int) len]; //crea array content del tamaño del canal de del fis
		int bytesRead = 0;
		int numRead = 0;
		String line;
		while ((line = fid.readLine()) != null) {
			// Put every character in upper case
			line = line.toUpperCase();
			numRead = line.length(); //Crea variable con longitud de linea del archivo de entrada q estoy leyendo. el numero de valores leidos será el tamaño de la linea que acaba de leer
			/*Creo array auxiliar que lleno con los datos (bytes de info) de las líneas del archivo de entrada. 
			 * en cada posición una secuencia/fila de nucleotidos del archivo de entrada
			 */
			byte[] newData = line.getBytes(); 
			for (int i = 0; i < numRead; i++) //ejecuto hasta el tamaño de la línea del archivo
				content[bytesRead + i] = newData[i]; //voy llenando el array content con los datos del array auxiliar. 
				//Escribo en content cada secuencia en la posición en la que ha acabado la anterior.			
				bytesRead += numRead; //va a indicarme la posición en la q ha acabado de escribir la ultima secuencia.
		}
		fis.close();
		this.content = content; //inicializa el array content
		this.validBytes = bytesRead; //inicializa el atributo validBytes
	}

	/**
	 * Provides the data array that contains the complete sequence of nucleotids
	 * extracted from the FASTA file.
	 * 
	 * @return The data array with each nucleotid taking one byte.
	 */
	public byte[] getContent() {
		return content;
	}

	/**
	 * Provides the amount of bytes in the data array that are valid. Since this
	 * array is created before the amount of bytes in the FASTA file that contain
	 * actual nucleotids are know, a worst-case scenario is assumed. So, only
	 * positions between content[0] and content[validBytes -1] have actual genomic
	 * data.
	 * 
	 * @return The number of valid bytes.
	 */
	public int getValidBytes() {
		return validBytes;
	}

	/**
	 * Returns the sequence of nucleotides of the provided size found at the
	 * provided position of the data array. If the initialPos + size is after the
	 * valid bytes of the array, it returns null.
	 * 
	 * @param initialPos The first character of the sequence.
	 * @param size       The length of the sequence.
	 * @return An String representing the sequence.
	 */
	public String getSequence(int initialPos, int size) {
		if (initialPos + size >= validBytes)
			return null;
		return new String(content, initialPos, size);
	}

	/*
	 * Helper method that checks if a pattern appears at a specific position in the
	 * data array. It _checks every byte of the pattern one by one_. If the pattern
	 * length would need to access a position after the valid bytes of the array, it
	 * throws a new FASTAException to indicate this fact.
	 * 
	 * Returns true if the pattern is present at the given position; false
	 * otherwise.
	 */
	private boolean compare(byte[] pattern, int position) throws FASTAException {
		if (position + pattern.length > validBytes) {
			throw new FASTAException("Pattern goes beyond the end of the file.");
		}
		boolean match = true;
		for (int i = 0; i < pattern.length; i++) { //recorre ENTERO el pattern mirando si cada elemento coincide con el de esa posición del array "content". 
			if (pattern[i] != content[position + i]) {
				match = false;
			}
		}
		return match; //la variable ha ido cambiando para cada elemento
	}

	/*
	 * Improved version of the compare method that stops checking elements of the
	 * pattern when one has been found to be different.
	 */
	private boolean compareImproved(byte[] pattern, int position) throws FASTAException {
	//HECHO POR MI	
		if (position + pattern.length > validBytes) {
			throw new FASTAException("Pattern goes beyond the end of the file.");
		}
		boolean match = true;
		for (int i = 0; i < pattern.length; i++) { //recorre el pattern mirando si cada elemento coincide con el de esa posición del array "content". 
			if (pattern[i] != content[position + i]) {
				match = false;
				break; //IMPROVEMENT: en cuanto una base del patrón no coincide, devuelve false
			}
		}
		return match; //la variable ha ido cambiando para cada elemento
	}

	/*
	 * Improved version of the compare method that returns the number of bytes in
	 * the pattern that are different from the ones present in the data array at the
	 * given position.
	 * 
	 * Returns the number of characters in the pattern that are different from the
	 * ones present in the indicated position.
	 */
	private int compareNumErrors(byte[] pattern, int position) throws FASTAException {
		//HECHO POR MI
		if (position + pattern.length > validBytes) {
			throw new FASTAException("Pattern goes beyond the end of the file.");
		}
		int diference = 0; //inicializo variable con la que voy a contar el numero de bases que difieren entre el pattern y el array de datos en la posicion dada.
		for(int i = 0; i < pattern.length; i++) { //recorro el pattern
			if (pattern[i] != content[position + i]) { //miro indice a indice si coincide con el array de datos 
				//Si difiere la base de ese indice, contabilizo que hay una diferencia
				diference++; //estoy contando todas las diferencias que haya. Quizá sería mejor si creo un break cuando diference>1, así deja de buscar. Para SNV solo contamos posicion si difiere como maximo en 1 base.
			}
		}
		return diference;
	}
	

	/**
	 * Implements a linear search to look for the provided pattern in the data
	 * array. Returns a List of Integers that point to the initial positions of all
	 * the occurrences of the pattern in the data.
	 * 
	 * @param pattern The pattern to be found.
	 * @return All the positions of the first character of every occurrence of the
	 *         pattern in the data.
	 */
	public List<Integer> search(byte[] pattern) {
	//HECHO POR MI
		List<Integer> initialPosition = new ArrayList<Integer>(); //Creo la lista que voy a devolver con las posiciones iniciales en las que aparece el patrón 
		for(int i = 0; i < this.content.length; i++) { //recorro el array de datos
			try {
				if (compareImproved (pattern, i)) { //ejecutando compare entre el pattern y cada posicion del array de datos
					initialPosition.add(i);
				}	
			} catch (FASTAException e) { //Tengo en cuenta que compare puede tirar FASTAException. Si se produce, la capturamos y terminamos la búsqueda.
				break;
			}
		}
		return initialPosition; //devuelvo la lista de las posiciones en las que he encontrado el pattern
	}

	/**
	 * Implements a linear search to look for the provided pattern in the data array
	 * but allowing a SNV (Single Nucleotide Variant). In SNV, one nucleotide is
	 * allowed to be different from the pattern. Therefore, this method returns a
	 * List of Integers that point to the initial positions of all the occurrences
	 * of the pattern in the data and all the occurrences of the pattern with one
	 * error in the data
	 * 
	 * @param pattern The pattern to be found.
	 * @return All the positions of the first character of every occurrence of the
	 *         pattern (with up to 1 errors) in the data.
	 */
	public List<Integer> searchSNV(byte[] pattern) {//Va a usar compareNumErrors
	//HECHO POR MI
		List<Integer> initialPosition = new ArrayList<Integer>(); //Creo la lista que voy a devolver con las posiciones iniciales en las que aparece el patrón permitiendo SNV
		for(int i = 0; i < this.content.length; i++) { //recorro el array de datos
			try {
				if (compareNumErrors(pattern, i) <= 1) { //si pattern difiere en 0 o en 1 con esa posicion del genoma, añado la posicion a mi lista de resultados.
					initialPosition.add(i);
				}	
			} catch (FASTAException e) { //Tengo en cuenta que compare puede tirar FASTAException. Si se produce, la capturamos y terminamos la búsqueda.
				break;
			}
		}
		return initialPosition; //devuelvo la lista de las posiciones en las que he encontrado el pattern permitiendo SNV
	}

	public static void main(String[] args) {
		long t1 = System.nanoTime();
		FASTAReader reader = new FASTAReader(args[0]);
		if (args.length == 1)
			return;
		System.out.println("Tiempo de apertura de fichero: " + (System.nanoTime() - t1));
		long t2 = System.nanoTime();
		List<Integer> posiciones = reader.searchSNV(args[1].getBytes());
		System.out.println("Tiempo de búsqueda: " + (System.nanoTime() - t2));
		if (posiciones.size() > 0) {
			for (Integer pos : posiciones)
				System.out.println("Encontrado " + args[1] + " en " + pos);
		} else
			System.out.println("No he encontrado : " + args[1] + " en ningun sitio");
		System.out.println("Tiempo total: " + (System.nanoTime() - t1));
	}
}

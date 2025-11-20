package es.upm.dit.aled.lab3.binary;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import es.upm.dit.aled.lab3.FASTAReader;

/**
 * Reads a FASTA file containing genetic information and allows for the search
 * of specific patterns within these data. The information is stored as an array
 * of bytes that contain nucleotides in the FASTA format. Since this array is
 * usually created before knowing how many characters in the origin FASTA file
 * are valid, an int indicating how many bytes of the array are valid is also
 * stored. All valid characters will be at the beginning of the array.
 * 
 * This extension of the FASTAReader uses a sorted dictionary of suffixes to
 * allow for the implementation of binary search.
 * 
 * @author mmiguel, rgarciacarmona
 *
 */
public class FASTAReaderSuffixes extends FASTAReader {
	protected Suffix[] suffixes;

	/**
	 * Creates a new FASTAReader from a FASTA file.
	 * 
	 * At the end of the constructor, the data is sorted through an array of
	 * suffixes.
	 * 
	 * @param fileName The name of the FASTA file.
	 */
	public FASTAReaderSuffixes(String fileName) {
		// Calls the parent constructor
		super(fileName);
		this.suffixes = new Suffix[validBytes];
		for (int i = 0; i < validBytes; i++)
			suffixes[i] = new Suffix(i);
		// Sorts the data
		sort();
	}

	/*
	 * Helper method that creates a array of integers that contains the positions of
	 * all suffixes, sorted alphabetically by the suffix.
	 */
	private void sort() {
		// Instantiate the external SuffixComparator, passing 'this' (the reader)
		// so it can access the content and validBytes fields.
		SuffixComparator suffixComparator = new SuffixComparator(this);
		// Use the external Comparator for sorting.
		Arrays.sort(this.suffixes, suffixComparator);
	}

	/**
	 * Prints a list of all the suffixes and their position in the data array.
	 */
	public void printSuffixes() {
		System.out.println("-------------------------------------------------------------------------");
		System.out.println("Index | Sequence");
		System.out.println("-------------------------------------------------------------------------");
		for (int i = 0; i < suffixes.length; i++) {
			int index = suffixes[i].suffixIndex;
			String ith = "\"" + new String(content, index, Math.min(50, validBytes - index)) + "\"";
			System.out.printf("  %3d | %s\n", index, ith);
		}
		System.out.println("-------------------------------------------------------------------------");
	}

	/**
	 * Implements a binary search to look for the provided pattern in the data
	 * array. Returns a List of Integers that point to the initial positions of all
	 * the occurrences of the pattern in the data.
	 * 
	 * @param pattern The pattern to be found.
	 * @return All the positions of the first character of every occurrence of the
	 *         pattern in the data.
	 */
	@Override
	public List<Integer> search(byte[] pattern) {
		// HECHO POR MI:
		List<Integer> positionFound = new ArrayList<Integer>(); //Creo la lista que voy a devolver con las posiciones en las que aparece el patrón
		// Inicializacion:
		int lo = 0; //Limite inferior de la busqueda binaria 
		int hi = suffixes.length; //Limite superior de la busqueda binaria
		boolean found = false; //Lo usaremos para determinar si se ha encontrado el pattern
		int index = 0; //contador que rastrea el carácter actual que se está comparando con el pattern, indice con el que recorremos el pattern y el sufijo dentro de content
		
		//Comparación iterativa (bucle de búsqueda binaria)
		do {
			int m = (int) Math.floor(lo+(hi-lo)/2); //Calculo el indice medio en el rango de busqueda actual
			int posSuffix = suffixes[m].suffixIndex; //Extraigo la posición en el genoma (suffixIndex) del sufijo que se encuentra en suffix[m]
			
			//Comienza a comparar el pattern con este sufijo carácter por carácter, comenzando por pattern[index]
			if(pattern[index] == content[posSuffix+index]) { //Si los caracteres actuales coinciden:
					
				index++; //Incremento index para verificar el siguiente caracter en la siguiente iteración
				//DUDA: incremento index antes o después: Tiene que ser ANTES, sino te sales del array.
				if(index == pattern.length){//COINCIDENCIA COMPLETA ENCONTRADA: comparación llega al final del patrón y los ultimos caracteres también coinciden.
					positionFound.add(posSuffix); //guardo la posición actual en la lista de resultados
					found = true;
					//MODIFICACIÓN PARA ENCONTRAR TODAS LAS COINCIDENCIAS: Cuando encuentro una, recorrer la lista de sufijos hacia arriba y hacia abajo.
					
					//Recorro lista de sufijos HACIA ARRIBA hasta que encuentro uno que no coincida
					int index2 = 0;//contador para recorrer patrón y los sufijos anteriores a la coincidencia hallada
					for (int i = 1; i<=m; i++) {
						int posSuffixUp = suffixes[m-i].suffixIndex;
						if(pattern[index2] != content[posSuffixUp]) { //si el primer valor del patrón no coincide con el primero del sufijo, dejo de mirar ese sufijo ni anteriores
							break;
						}
						while(pattern[index2] == content[posSuffixUp+index2]) {
							index2++;
							if(index2 == pattern.length) {
								positionFound.add(posSuffixUp);
								index2=0;
								break;
							}
						}
					}

					//Recorro lista de sufijos HACIA ABAJO hasta que encuentro uno que no coincida
					int index3 = 0;//contador para recorrer patrón y los sufijos posteriores a la coincidencia hallada
					for (int i = 1; i<(suffixes.length -m) ; i++) {
						int posSuffixDw = suffixes[m+i].suffixIndex;
						if(pattern[index3] != content[posSuffixDw]) { //si el primer valor del patrón no coincide con el primero del sufijo, dejo de mirar ese sufijo ni anteriores
							break;
						}
						while(pattern[index3] == content[posSuffixDw+index3]) {
							index3++;
							if(index3 == pattern.length) {
								positionFound.add(posSuffixDw);
								index3=0;
								break;
							}
						}
					}	
				}
			}
			//DIVISION ESTANDAR DE BUSQUEDA BINARIA: si el principio del sufijo de suffixes[m] no coincide con el patrón
			else if (pattern[index] < content[posSuffix+index]) { // el caracter del pattern es lexicográficamente anterior al carácter del sufijo
				hi = m--; //restablecemos el limite superior de la busqueda descartando la mitad derecha de suffixes
				index = 0; //IMP: reinicio el contador
			}
			else {//pattern[index] > content[posSuffix+index]: el caracter del pattern es lexicográficamente posterior al del sufijo
				lo = m++; //restablecemos el limite inferior de la busqueda descartando la mitad izquierda de suffixes
				index = 0; //IMP: reinicio el contador
			}
			//RESTABLECEMOS LIMITES SUPERIOR E INFERIOR sin coger de nuevo m, que ya lo hemos mirado.
		} while (!found==true && !(hi-lo<=1)); //Para que el bucle continúe hasta que se haya encontrado una coincidencia total O el espacio de busqueda sea demasiado pequeño
			
		return positionFound;//devuelvo la lista de las posiciones en las que he encontrado el pattern		
	}

	
	public static void main(String[] args) {
		long t1 = System.nanoTime();
		FASTAReaderSuffixes reader = new FASTAReaderSuffixes(args[0]);
		if (args.length == 1)
			return;
		byte[] patron = args[1].getBytes();
		System.out.println("Tiempo de apertura de fichero: " + (System.nanoTime() - t1));
		long t2 = System.nanoTime();
		System.out.println("Tiempo de ordenación: " + (System.nanoTime() - t2));
		reader.printSuffixes();
		long t3 = System.nanoTime();
		List<Integer> posiciones = reader.search(patron);
		System.out.println("Tiempo de búsqueda: " + (System.nanoTime() - t3));
		if (posiciones.size() > 0) {
			for (Integer pos : posiciones)
				System.out.println("Encontrado " + args[1] + " en " + pos);
		} else
			System.out.println("No he encontrado " + args[1] + " en ningún sitio.");
		System.out.println("Tiempo total: " + (System.nanoTime() - t1));
	}
}

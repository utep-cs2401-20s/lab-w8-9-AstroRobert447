class AminoAcidLL{
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;
  


  AminoAcidLL(){

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon 
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon){
	  aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
	  codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
	  counts = new int[codons.length];
	  inerCodon(inCodon);
	  next = null;
  }
  
  private void inerCodon(String inCodon) {
	  for(int i = 0; i<codons.length; i++) {
		  if(codons[i].equals(inCodon.toUpperCase())) {
			  counts[i]++;
		  }
		  
		  else {
			  counts[i] = 0;
		  }
	  }
  }

  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
  private void addCodon(String inCodon){
	  //What you see here is here is the code that is relate with the one that you showing us in the lab class.
	  //However I keep getting error because somehow it keep asking me to add the stuff in the imports.
	  
	  //Here is the example of the original.
	  
	  
	 // if(next == null) {
	//	  if(node == codons) {
	//		  int increCodons(AminoAcidResources.getAminoAcidFromCodon(codons = aminoAcid));
	//	  }
	//	  if(aminoAcid == AminoAcidResources.getCodonListForAminoAcid(c)) {
	//		  inerCodons(c);
	//	  }
	// }
	  
	  //Therefore I decided to figure out in my own way and this is what I made. The instructor said its was okay
	  //doing in my own way.
	  char c = AminoAcidResources.getAminoAcidFromCodon(inCodon);
	  String[] str =  AminoAcidResources.getCodonListForAminoAcid(c);
	  for(int i = 0; i<str.length; i++) {
		  if(inCodon.equals(str[i])) {
			  this.inerCodon(inCodon);
		  }
		  else if(this.next != null) {
			  this.next.addCodon(inCodon);
		  }
		  else {
			  this.next = new AminoAcidLL(inCodon);
		  }
	  
	  }
	  
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;
    for(int i = 0; i < counts.length; i++) {
    	sum += counts[i];
    }
    return sum;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    return codonDiff(inList);
  }

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList){
    return totalDiff(inList);
  }


  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList(){
    char[] check = null;
    for(int i = 0; i < codons.length; i++) {
    	check [i] = AminoAcidResources.getAminoAcidFromCodon(codons[i]);
    }
    return check;
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    return counts;
  }


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    if (this.next == null) {
    	return true;
    }
    if(this.aminoAcid > this.next.aminoAcid) {
    	return false;
    }
    else if(this.next.next != null) {
    	this.next.isSorted();
    }
    return true;
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
	  String firstThree = inSequence.substring(0,3);
	  AminoAcidLL firstNode = new AminoAcidLL(firstThree);
	  char AA = AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0,3));
	  if(AA == '*') {
		return firstNode;
	  }
	  for(int i = 3; i < inSequence.length(); i += 3) {
		  firstNode.addCodon(inSequence.substring(i, i+3));
		  AA = AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(i, i+3));
		  if(AA == '*') {
			  break;
		  
		  }
	  }
	  return firstNode;
  }


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){
    if(inList == null || inList.next == null) {
    	return inList;
    }
    
    AminoAcidLL lowSpeed = inList, highSpeed = inList;
    	while (highSpeed.next != null && highSpeed.next.next != null) {//Right here, the instructor had
    		lowSpeed = lowSpeed.next;								   //explain the different between the While and For loops
    		highSpeed = highSpeed.next.next;						   //and understand why the For did not want to work on the loop.
    	}
    	
    	AminoAcidLL nextToMiddlePointer = lowSpeed.next;
    	lowSpeed.next = null;
    	
    	AminoAcidLL First = sort(inList);
    	AminoAcidLL Second = sort(nextToMiddlePointer);
    	AminoAcidLL sortedList = merge(First, Second);
    	return sortedList;
    	
  }
  
  private static AminoAcidLL merge(AminoAcidLL First, AminoAcidLL Second) {
	  if (First == null) {
		  return Second;
	  }
	  if (Second == null) {
		  return First;
	  }
	  
	  AminoAcidLL result = null;
	  
	  if(First.aminoAcid <= Second.aminoAcid) {
		  result = First;
		  result.next = merge(First.next, Second);
	  }
	  else {
		  result = First;
		  result.next = merge(First, Second.next);
	  }
	  return result;
  }
}
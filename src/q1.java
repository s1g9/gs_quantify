import java.util.*;
import java.math.*;


class q1 {

	public static void main(String[] args) {
		Scanner scan = new Scanner(System.in);
		int n;
		String firstline = new String();
		firstline = scan.nextLine();
		n = Integer.parseInt(firstline.split(" ")[1]);
		String line= new String();
		HashMap<Integer, HashMap<String, ArrayList<Integer>>> outermap = new HashMap<Integer, HashMap<String, ArrayList<Integer>>>();
		HashMap<String, ArrayList<Integer>> innermap = new HashMap<String, ArrayList<Integer>>();
		int previous=-1;
		ArrayList<Integer> timestamps=new ArrayList<Integer>();
		for(int i=0; i<n; i++) {
		line = scan.nextLine();
		if(Integer.parseInt(line.split(" ")[0])!=previous){
		timestamps.add(Integer.parseInt(line.split(" ")[0]));
	    innermap = new HashMap<String, ArrayList<Integer>>();
		}
		
		for(int j=2; j < line.split(" ").length; j=j+2){
			ArrayList<Integer> fieldvalues=new ArrayList<Integer>();
			if(innermap.get(line.split(" ")[1].concat("_"+line.split(" ")[j]))!=null){fieldvalues=innermap.get(line.split(" ")[1].concat("_"+line.split(" ")[j]));}
			fieldvalues.add(Integer.parseInt(line.split(" ")[j+1]));
			innermap.put(line.split(" ")[1].concat("_"+line.split(" ")[j]), fieldvalues);
		}
		outermap.put(Integer.parseInt(line.split(" ")[0]), innermap);
		previous=Integer.parseInt(line.split(" ")[0]);
		}
		System.out.println("tickfile completed");
		
		while(scan.hasNextLine()){
			line=scan.nextLine();
			switch (line.split(" ")[0]) {
            case "sum": summation(timestamps, outermap,Integer.parseInt(line.split(" ")[1]),Integer.parseInt(line.split(" ")[2]),line.split(" ")[3],line.split(" ")[4]);
                     break;
            case "product": productof(timestamps, outermap,Integer.parseInt(line.split(" ")[1]),Integer.parseInt(line.split(" ")[2]),line.split(" ")[3],line.split(" ")[4],line.split(" ")[5] ) ;
                     break;
            case "max": maxof(timestamps, outermap,Integer.parseInt(line.split(" ")[1]),Integer.parseInt(line.split(" ")[2]),line.split(" ")[3],line.split(" ")[4],Integer.parseInt(line.split(" ")[5]) ) ;
                     break;
            case "delta": deltaof(timestamps, outermap,line.split(" ")[1],line.split(" ")[2],Integer.parseInt(line.split(" ")[3])) ;
                     break;
            default:  System.out.println("Invalid function asked");
                     break;
        }
		}
	}


	private static void deltaof(
			ArrayList<Integer> timestamps, HashMap<Integer, HashMap<String, ArrayList<Integer>>> outermap, String symbol, String field, int k) {
		ArrayList<Integer> time=new ArrayList<Integer>();
		ArrayList<Integer> valarray=new ArrayList<Integer>();
		for (int i = 0; i<timestamps.size(); i++){
			Map<String, ArrayList<Integer>> outmap = outermap.get(timestamps.get(i));
			ArrayList<Integer> value = outmap.get(symbol.concat("_"+field));
			if(value==null){continue;} else {
					valarray.add(value.get(0));
					time.add(timestamps.get(i));
		}
	}
	System.out.println(leastsegmemtationerror(time,valarray,k));
}




	private static int leastsegmemtationerror(ArrayList<Integer> timestamps,
			ArrayList<Integer> valarray, int k) {
		int no_of_segments=(int) Math.ceil(((double)timestamps.size())/2);
		int errorbottom,errortop,error,minerror;
		errorbottom=k*no_of_segments;
		int leasterror=(int) square_error(timestamps,valarray);
		errortop=k+leasterror+1;
		minerror=errortop;
		ArrayList<Integer> err = new ArrayList<Integer>();
		ArrayList<Integer> a = new ArrayList<Integer>();
		ArrayList<Integer> b= new ArrayList<Integer>();
		for(int i=0; i <k;i=i+2){
			a=(ArrayList<Integer>) timestamps.subList(i,timestamps.size()-1);
			b=(ArrayList<Integer>) valarray.subList(i,timestamps.size()-1);
			error=(int) square_error(a,b);
			err.add((k-i)+error+1);
			minerror=Math.min(minerror,(k-i)+error+1);
		}
		return Math.min(errorbottom,errortop);
	}




	private static void maxof(ArrayList<Integer> timestamps, HashMap<Integer, HashMap<String, ArrayList<Integer>>> outermap, int starttime,
			int endtime, String symbol, String field, int k) {
		ArrayList<Integer> valarray=new ArrayList<Integer>();
		int end = binarySearch(timestamps, endtime,0);
		for (int i = binarySearch(timestamps, starttime,1); i<=end; i++){
			Map<String, ArrayList<Integer>> outmap = outermap.get(timestamps.get(i));
			ArrayList<Integer> value = outmap.get(symbol.concat("_"+field));
			if(value==null){continue;} else {valarray.addAll(value);}
		}
		Collections.sort(valarray);
		for(int i=0;i<k && i< valarray.size(); i++){
			System.out.print(valarray.get(valarray.size()-1-i) +" ");	
		}
		System.out.println();
	}


	private static void productof(ArrayList<Integer> timestamps, HashMap<Integer, HashMap<String, ArrayList<Integer>>> outermap, int starttime,
			int endtime, String symbol, String field1, String field2) {
		int sum=0;
		int end =  binarySearch(timestamps, endtime,0);
		for (int i = binarySearch(timestamps, starttime,1); i<=end; i++){
			Map<String, ArrayList<Integer>> outmap = outermap.get(timestamps.get(i));
			ArrayList<Integer> value1 = new ArrayList<Integer>();
			if(outmap.get(symbol.concat("_"+field1))!=null){
			value1.addAll( outmap.get(symbol.concat("_"+field1)));}
			ArrayList<Integer> value2 = new ArrayList<Integer>();
			if(outmap.get(symbol.concat("_"+field2))!=null){
			value2.addAll( outmap.get(symbol.concat("_"+field2)));}
			int value=0;
			if(value1.isEmpty() || value2.isEmpty()){
				continue;
			}
			for(int j=0; j< Math.min(value1.size(), value2.size()) ; j++) {
				value+=value1.get(j)*value2.get(j);
			}
			sum=sum+value;
		}
		System.out.println(sum);
		
	}


	private static void summation(
			ArrayList<Integer> timestamps, HashMap<Integer, HashMap<String, ArrayList<Integer>>> outermap, int starttime,
			int endtime, String symbol, String field) {
		int sum=0;
		int end =  binarySearch(timestamps, endtime,0);
		for (int i = binarySearch(timestamps, starttime,1); i<=end; i++){
			Map<String, ArrayList<Integer>> outmap = outermap.get(timestamps.get(i));
			ArrayList<Integer> value = new ArrayList<Integer>();
			if(outmap.get(symbol.concat("_"+field))!=null){
			value.addAll( outmap.get(symbol.concat("_"+field)));}
			int sumvalue=0;
			if(value.isEmpty()){
				continue;
			}
			for(int j=0; j< value.size() ; j++) {
				sumvalue+=value.get(j);
			}
			sum = sum + sumvalue ;
		}
		System.out.println(sum);
	}

	private static int binarySearch(ArrayList<Integer> data, int key, int upper) 
	    {
	         int low = 0;
	         int high = data.size() - 1;
	          
	         while(high >= low) {
	             int middle = (low + high) / 2;
	             if(data.get(middle) == key) {
	                 return middle;
	             }else if(data.get(middle) < key) {
	                 low = middle + 1;
	             }else if(data.get(middle) > key) {
	                 high = middle - 1;
	             }
	        }
	         if(upper==1){
	        	 return low;
	         } else return high;
	    }
	public static double variance(ArrayList<Integer> v) {
	    double mu = mean(v);
	    double sumsq = 0.0;
	    for (int i = 0; i <v.size(); i++)
	      sumsq += Math.pow(mu - v.get(i),2);
	    return sumsq / (v.size());

	  }
	public static double mean(ArrayList<Integer> v) {
	    double tot = 0.0;
	    for (int i = 0; i <v.size(); i++)
	      tot += v.get(i);
	    return tot / (v.size());
	  }
	public static double covar(ArrayList<Integer> v1, ArrayList<Integer> v2) {
	    double m1 = mean(v1);
	    double m2 = mean(v2);
	    double sumsq = 0.0;
	    for (int i = 0; i <v1.size(); i++)
	      sumsq += (m1 - v1.get(i)) * (m2 - v2.get(i));
	    return sumsq / (v1.size());
	  }
	public static double square_error(ArrayList<Integer> v1, ArrayList<Integer> v2) {
		double a,b,error=0;
		double err = 0; 
		a=covar(v1,v2)/variance(v1) ;
		b=mean(v2)-a*mean(v1);
//		err.add(0);// Initializing the total error i.e first element
//		err.add(0);// Initializing the Maximum error i.e. second element
//		ArrayList<Integer> maxindexes = new ArrayList<Integer>(); 
		for(int j=0;j<v1.size();j++){
			error=Math.pow((v2.get(j)-a*v1.get(j)-b),2);
//			if(error>err.get(1)){
//				ArrayList<Integer> maxindexes = new ArrayList<Integer>(); 
//				maxindexes.add(i);
//			}
//			err.get(1)=Math.max(err.get(1),error);
//			if(err.get(1)==err && maxindex.get(0)!=i){
//			maxindexs.add(i);// Storing indices of Maximum error points}
			err=err+error;
		}
//		err.addall(maxindexes);
	    return err; // Returning the total error, Maximum error and array of indices corresponding to maximum error respectively
	  }
}



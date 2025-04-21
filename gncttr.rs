#![allow(unstable)] // allow unstable libraries
#![feature(rustc_private)] 

#[macro_use] extern crate lazy_static;

extern crate regex;

use std::io;
use std::io::prelude::*;
//use chomp::primitives::Primitives;
use std::collections::HashMap;
use regex::Regex;
use regex::Error;


/// write: how compiled? how run programm? which parameters for programm call?
// cat gncttr/*fna | cargo run --

// programm for extraction of codon sequences from genoms
fn main () {
	// fill hash with genetic code
	let mut code = HashMap::new();
	code = gencode(); // HashMap<String, String>
	let starts : [String; 3] = [format!("ATG"), format!("TTG"), format!("CTG")]; // array -> fixed size!
	let limit = 100000; // length limitation for analysed sequence window
	
	// regex for identication of nucleodide
//	let nc_pattern = r"[ATCGN]";
//	let res = Regex::new(nc_pattern);
	
	// translate codon into AS
/*	let t = "TTT".to_string();
	let AS = code.get(&t);
	println!("AS for {} is {}", t, AS.unwrap());
	if let Ok(regex) = res {
		println!("regex ok, trying to match..");
		let matches = regex.is_match(&t);
		if matches {
			println!("{} is nucleotid {}", t, nc_pattern);
		} else {
			println!("{} NOT a nucleotid sequence {}", t, nc_pattern);
		}
	} else {
		println!("{:?}", res);
	}
*/	
	
	// playing with regex
	let nc_pattern = Regex::new(r"[ATCGN]").unwrap(); // assign regex
	let head = Regex::new(r"^>(.*)").unwrap(); // regex for capturing header
	let test = ">wppwpweATTCGTGC8023582305";
//	let after = re.replace_all(test, "--");
//	println!("{}", after);
/*	for cap in re.captures_iter(test) {
		println!("{} -> {} -> {}", &cap[0], &cap[1], &cap[2]); // &cap[0] -> all captured by regex, &cap[1] -> 1st captured group, &cap[2] -> 2nd captured group, ...
	}
*/
/*	for cap in head.captures_iter(test) {
		println!("Header: {}", &cap[1]);
	}
*/

	//let matches = re.is_match(&test);
	//assert!(re.is_match("TCGN"));
	
	// read sequence from stdin
	let mut input = io::stdin();
	let mut linebuffer = String::new();
	loop {
		match input.read_line(&mut linebuffer) {
			Ok(i) if i == 0 => { break; },
			Ok(i) => {
				//println!("{} {}", linebuffer, i);
				let h = head.is_match(&linebuffer); // fasta header?
				let nt = nc_pattern.is_match(&linebuffer); // nucleotide sequences?
				if h {
					//eprint!("HEADER found! {}", linebuffer);
				} else if nt {
					eprint!("NUCLEOTIDE sequences {}", linebuffer);
					
					// split reads into triples -> remains shift into next reads
			
					// startcodon?	
					//println!("{}", i%3);
					let res = i%3;
					println!("last 2: {:?}", slice_from_end(&linebuffer, 2));
					println!("first 2: {:?}", slice_from_start(&linebuffer, 2));
					// test if codon is startcodon -> TRUE: save sequence
				}
			},
			Err(e) => { panic!(e); }
		}
		linebuffer.clear();
	}

	
	//let stdin = io::stdin();
	//let mut seq = String::new();
	//match read_line(&mut stdin) {
		//Ok(n) => {
			//println!("{} bytes read", n);
			//println!("{}", read_line(&stdin));
		//}
		//Err(error) => println!("error: {}", error);
	//}

/*	for mut line in stdin.lock().lines() {
		let mut read = line.unwrap().to_string();
		let h = head.is_match(&read); // fasta header?
		
		if h {
			//eprintln!("HEADER found! {}", read);
		} else if nt {
			//eprintln!("NUCLEOTIDE sequences {}", read);
			let res = read.len() % 3;
			println!("overlap: {}", res);
			
			// split reads into triples -> remains shift into next reads
			
			// test if codon is startcodon -> TRUE: save sequence
			
			chomp(&mut read);
			println!("{}", read.len()); // 80 nt
				//seq.push_str(&read);
/*				if seq.len() > limit {
					
					let duple = seq.clone();
					let codons = triple(&seq);
					//println!("{}", seq.len());
					seq = seq[81..(seq.len() - 81)].to_string();
				}
*/				
				// search for startcodon in seq
				
				// search for stopcodon in seq
				
				// start & stopcodon?? -> get sequence!

			//seq.push_str(&read); // concatenate strings
			//println!("{}", seq);
		} else {
			eprintln!("[ERROR]\tskipping line because of non-nucleotid sequence!");
			continue;
		}
	}
*/

	// split sequence into codons
/*	let codons = triple(&seq);
	for &c in codons.iter() {
		//println!("{}", c);
		
		// find startcodon
		//let ini = starts.iter().position(|r| r == c).unwrap();
		//println!("{}", ini);
		
		// translate codons
		let AS = code.get(c);
		match AS {
			Some(AS) => println!("{}", AS),
			None 	=> println!("Something went wrong!!\n{} is no codon!", c),		  //println!("AS for codon {} is {}!", c, AS);
		}

	}
*/	
	
	// clean sequence -> remove header (![ATCGNU])

/*	// remove newline
	let mut test = format!("test\n");
	let dummy = format!("dummy1");
	//println!("{}", test);
	//rm_newline(&test);
	//test.pop();
	chomp(&mut test);
	//println!("{}", test);
	test.push_str(&dummy);
	println!("{}", test);
*/
	// split sequence into triplet
	
	// compare & translate triplet into AS
	
	// CDS??
	
	
	// print out CDS


}

// remove n first character of string
fn slice_from_start(s: &str, n: usize) -> Option<&str> {
    s.char_indices().nth(n).map(|(i, _)| &s[i..])
}

// slice n last character from a string
fn slice_from_end(s: &str, n: usize) -> Option<&str> {
    s.char_indices().rev().nth(n).map(|(i, _)| &s[i..])
}

// remove newline character from string
fn triple (string: &str) -> Vec<&str> {
	// split into triplets
	let mut triplet = String::new();
	let mut codon: Vec<&str> = Vec::new();;
	for (i, s) in string.chars().enumerate() {
		if i%3 == 0 {
			codon.push(&string[(i)..(i+3)]);
			triplet = format!(""); // empty string
		}
		triplet.push(s);
	}
	codon
}


// remove newline character at end of string // rebuild of trim_right_matches("\n")
fn chomp(raw: &mut String) {
	let last = raw.chars().last().unwrap();
	//println!("last character of string is: {}", last);
	if last == '\n' {
		eprintln!("[INFO]\tNewline removed!");
		raw.pop();
	};
}


// create hash containing genetic code
fn gencode() -> (HashMap<String, String>) {
	let mut code = HashMap::new();
	code.insert("TTT".to_string(), "F".to_string());
	code.insert("TTC".to_string(), "F".to_string());
	code.insert("TTA".to_string(), "L".to_string()); 
	code.insert("TTG".to_string(), "L".to_string()); 
	code.insert("TCT".to_string(), "S".to_string());
	code.insert("TCC".to_string(), "S".to_string());
	code.insert("TCA".to_string(), "S".to_string());
	code.insert("TCG".to_string(), "S".to_string());
	code.insert("TAT".to_string(), "Y".to_string());
	code.insert("TAC".to_string(), "Y".to_string());
	code.insert("TAA".to_string(), "*".to_string());
	code.insert("TAG".to_string(), "*".to_string());
	code.insert("TGT".to_string(), "C".to_string());
	code.insert("TGC".to_string(), "C".to_string());
	code.insert("TGA".to_string(), "*".to_string());
	code.insert("TGG".to_string(), "W".to_string());
	code.insert("CTT".to_string(), "L".to_string());
	code.insert("CTC".to_string(), "L".to_string());
	code.insert("CTA".to_string(), "L".to_string());
	code.insert("CTG".to_string(), "L".to_string());
	code.insert("CCT".to_string(), "P".to_string());
	code.insert("CCC".to_string(), "P".to_string());
	code.insert("CCA".to_string(), "P".to_string());
	code.insert("CCG".to_string(), "P".to_string());
	code.insert("CAT".to_string(), "H".to_string());
	code.insert("CAC".to_string(), "H".to_string());
	code.insert("CAA".to_string(), "Q".to_string());
	code.insert("CAG".to_string(), "Q".to_string());
	code.insert("CGT".to_string(), "R".to_string());
	code.insert("CGC".to_string(), "R".to_string());
	code.insert("CGA".to_string(), "R".to_string());
	code.insert("CGG".to_string(), "R".to_string());
	code.insert("ATT".to_string(), "I".to_string());
	code.insert("ATC".to_string(), "I".to_string());
	code.insert("ATA".to_string(), "I".to_string());
	code.insert("ACT".to_string(), "T".to_string());
	code.insert("ACC".to_string(), "T".to_string());
	code.insert("ACA".to_string(), "T".to_string());
	code.insert("ACG".to_string(), "T".to_string());
	code.insert("AAT".to_string(), "N".to_string());
	code.insert("AAC".to_string(), "N".to_string());
	code.insert("ATG".to_string(), "M".to_string());
	code.insert("AAA".to_string(), "K".to_string());
	code.insert("AAG".to_string(), "K".to_string());
	code.insert("AGT".to_string(), "S".to_string());
	code.insert("AGC".to_string(), "S".to_string());
	code.insert("AGA".to_string(), "R".to_string());
	code.insert("AGG".to_string(), "R".to_string());
	code.insert("GTT".to_string(), "V".to_string());
	code.insert("GTC".to_string(), "V".to_string());
	code.insert("GTA".to_string(), "V".to_string());
	code.insert("GTG".to_string(), "V".to_string());
	code.insert("GCT".to_string(), "A".to_string());
	code.insert("GCC".to_string(), "A".to_string());
	code.insert("GCA".to_string(), "A".to_string());
	code.insert("GCG".to_string(), "A".to_string());
	code.insert("GAT".to_string(), "D".to_string());
	code.insert("GAC".to_string(), "D".to_string());
	code.insert("GAA".to_string(), "E".to_string());
	code.insert("GAG".to_string(), "E".to_string());
	code.insert("GGT".to_string(), "G".to_string());
	code.insert("GGC".to_string(), "G".to_string());
	code.insert("GGA".to_string(), "G".to_string());
	code.insert("GGG".to_string(), "G".to_string());
	return code;
}

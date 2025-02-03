use std::env;
use std::fs::File;
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;
use rayon::{self, ThreadPoolBuilder, prelude::*};

fn get_sequences_from_fasta(file: &String) -> Vec<String> {

    let file = File::open(file).expect("Unable to open file");
    let reader = BufReader::new(file);

    let mut sequences: Vec<String> = Vec::new();
    let mut current_sequence: String = String::new();

    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        if line.starts_with('>') {
            if !current_sequence.is_empty() {
                sequences.push(current_sequence.clone());
                current_sequence.clear();
            }
        } else {
            current_sequence.push_str(&line);
        }
    }
    if !current_sequence.is_empty() {
        sequences.push(current_sequence);
    }
    return sequences;
}

fn get_sequences_from_fastq(file: &String) -> Vec<String> {

    let file = File::open(file).expect("Unable to open file");
    let reader = BufReader::new(file);
    let sequences: Vec<String> = Vec::new();

    let mut current_sequence: String = String::new();

    let mut read: bool = false;

    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        if line.starts_with("@") {
            read = true;
        }
        if line.starts_with("+") {
            read = false;
            if !current_sequence.is_empty() {
                current_sequence.clear();
            }
        }
        if read {
            current_sequence.push_str(&line); 
        }
    }
    return sequences;
}

fn get_sequences(file: &String) -> Vec<String> {
    if file.ends_with("a") {
        let sequences: Vec< String> = get_sequences_from_fasta(file);
        return sequences;
    }
    if file.ends_with("q") {
        let sequences: Vec< String> = get_sequences_from_fastq(file);
        return sequences;
    }
    panic!("File must end with a or q to be recognized as FASTA or FASTQ");     
}

fn count_kmers(file: &String, k: usize) -> HashMap<String, u32> {

    // Shared HashMap protected by a Mutex for thread-safe updates
    let kmer_hashmap = Mutex::new(HashMap::new());

    // Get sequences from the file
    let sequences: Vec<String> = get_sequences(file);

    println!("Read {} sequences", sequences.len());
    println!("-------------------------------------");
    println!("Processed sequences:");

    // Atomic counter to track total processed sequences
    let progress = AtomicUsize::new(0);

    // Process sequences in parallel
    sequences.par_iter().enumerate().for_each(|(_i, sequence)| {

        // Local HashMap for each thread to reduce contention
        let mut local_map: HashMap<String, u32> = HashMap::new();

        // Increment the processed sequences counter
        let total_progress = progress.fetch_add(1, Ordering::Relaxed) + 1;

        if total_progress % (sequences.len()/10) == 0 {
            println!("{}", total_progress);
        }

        for j in 0..(sequence.len() - k) {
            *local_map.entry(sequence[j..j + k].to_string()).or_insert(0) += 1;
        }

        // Merge local HashMap into the global one
        let mut global_map = kmer_hashmap.lock().unwrap();
        for (key, value) in local_map {
            *global_map.entry(key).or_insert(0) += value;
        }
    });

    // Return the final HashMap
    Mutex::into_inner(kmer_hashmap).unwrap()
}

fn save_kmers(kmer_hashmap: HashMap<String, u32>) -> io::Result<()> {

    let mut file = File::create("kmer_counts.tsv")?;

    for (key, value) in &kmer_hashmap {
        writeln!(file, "{}\t{}", key, value)?;
    }

    Ok(())
}

fn main() -> io::Result<()> {

    let start = Instant::now();

    let args: Vec<String> = env::args().collect();

    let file: String = String::from(&args[1]);
    let k: usize  = args[2].parse::<usize>().unwrap();
    let threads   = args[3].parse::<usize>().unwrap();

    println!("-------------------------------------");
    println!("Arguments:");
    println!("File:    {}", file);
    println!("k:       {}", k);
    println!("Threads: {}", threads);
    println!("-------------------------------------");

    // Setup for parallel kmer counting
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    // Kmer counting
    let kmer_hashmap: HashMap<String, u32> = count_kmers(&file, k);

    println!("-------------------------------------");
    println!("Writing kmer counts to file");
    println!("-------------------------------------");

    save_kmers(kmer_hashmap).expect("File writing failed");

    let end = Instant::now();

    println!("DONE after {:?}", end.duration_since(start));

    Ok(())

}
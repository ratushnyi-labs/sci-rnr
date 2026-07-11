//! rnr -- production CLI for the RNR1 Type-I coder.
//!
//! Statically compiled port of the normative reference `impl/rnr1.py`
//! (archive format: impl/FORMAT.md).  The reference stays normative:
//! this binary is required to produce archives BYTE-IDENTICAL to
//! `rnr1.pack()` and to decode exactly the archives the reference
//! accepts (verified by rnr-cli/run_checks.py against the corpus
//! matrix).
//!
//! Subcommands (flags mirror the reference CLI, plus --threads):
//!
//!   rnr pack   <input> <output.rnr> [--W 3] [--K 65536] [--threads 1]
//!   rnr unpack <archive.rnr> <output> [--no-verify] [--threads 1]
//!   rnr read   <archive.rnr> --pos P --len k [--out FILE] [--verify]
//!   rnr info   <archive.rnr>

mod coder;
mod container;

use std::io::Write;
use std::process::ExitCode;

use container::{Archive, VERSION};

const DEFAULT_W: u32 = 3;
const DEFAULT_K: u32 = 65536;

fn die(msg: &str) -> ExitCode {
    eprintln!("rnr: error: {}", msg);
    ExitCode::from(1)
}

fn usage() -> ExitCode {
    eprintln!(
        "usage: rnr pack <input> <output> [--W {}] [--K {}] [--threads 1]\n\
         \x20      rnr unpack <archive> <output> [--no-verify] [--threads 1]\n\
         \x20      rnr read <archive> --pos P --len k [--out FILE] [--verify]\n\
         \x20      rnr info <archive>",
        DEFAULT_W, DEFAULT_K
    );
    ExitCode::from(2)
}

struct Args {
    positional: Vec<String>,
    w: u32,
    k: u32,
    threads: u32,
    no_verify: bool,
    verify: bool,
    pos: Option<u64>,
    len: Option<u64>,
    out: Option<String>,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut a = Args {
        positional: Vec::new(),
        w: DEFAULT_W,
        k: DEFAULT_K,
        threads: 1,
        no_verify: false,
        verify: false,
        pos: None,
        len: None,
        out: None,
    };
    let mut it = argv.iter();
    while let Some(arg) = it.next() {
        let mut take = |name: &str| -> Result<String, String> {
            it.next()
                .cloned()
                .ok_or_else(|| format!("missing value for {}", name))
        };
        match arg.as_str() {
            "--W" => a.w = take("--W")?.parse().map_err(|_| "bad --W")?,
            "--K" => a.k = take("--K")?.parse().map_err(|_| "bad --K")?,
            "--threads" => {
                a.threads = take("--threads")?.parse().map_err(|_| "bad --threads")?
            }
            "--no-verify" => a.no_verify = true,
            "--verify" => a.verify = true,
            "--pos" => a.pos = Some(take("--pos")?.parse().map_err(|_| "bad --pos")?),
            "--len" => a.len = Some(take("--len")?.parse().map_err(|_| "bad --len")?),
            "--out" => a.out = Some(take("--out")?),
            s if s.starts_with("--") => return Err(format!("unknown flag {}", s)),
            s => a.positional.push(s.to_string()),
        }
    }
    Ok(a)
}

fn read_file(path: &str) -> Result<Vec<u8>, String> {
    std::fs::read(path).map_err(|e| format!("{}: {}", path, e))
}

fn write_file(path: &str, data: &[u8]) -> Result<(), String> {
    std::fs::write(path, data).map_err(|e| format!("{}: {}", path, e))
}

fn cmd_pack(a: &Args) -> Result<(), String> {
    if a.positional.len() != 2 {
        return Err("pack needs <input> <output>".into());
    }
    let data = read_file(&a.positional[0])?;
    let raw = container::pack(&data, a.w, a.k, a.threads)?;
    write_file(&a.positional[1], &raw)?;
    let m = if data.is_empty() {
        0
    } else {
        (data.len() as u64 + a.k as u64 - 1) / a.k as u64
    };
    println!(
        "packed {} bytes -> {} bytes ({:.4} bpb), {} sub-block(s)",
        data.len(),
        raw.len(),
        8.0 * raw.len() as f64 / (data.len().max(1)) as f64,
        m
    );
    Ok(())
}

fn cmd_unpack(a: &Args) -> Result<(), String> {
    if a.positional.len() != 2 {
        return Err("unpack needs <archive> <output>".into());
    }
    let raw = read_file(&a.positional[0])?;
    let arc = Archive::parse(&raw)?;
    let data = arc.unpack(!a.no_verify, a.threads)?;
    write_file(&a.positional[1], &data)?;
    println!(
        "unpacked {} bytes -> {} bytes (verified: {})",
        raw.len(),
        data.len(),
        if a.no_verify { "no" } else { "yes" }
    );
    Ok(())
}

fn cmd_read(a: &Args) -> Result<(), String> {
    if a.positional.len() != 1 {
        return Err("read needs <archive> --pos P --len k".into());
    }
    let (pos, len) = match (a.pos, a.len) {
        (Some(p), Some(l)) => (p, l),
        _ => return Err("read needs --pos and --len".into()),
    };
    let raw = read_file(&a.positional[0])?;
    let arc = Archive::parse(&raw)?;
    let (data, warmup, decoded) = arc.read(pos, len, a.verify)?;
    match &a.out {
        Some(path) => write_file(path, &data)?,
        None => {
            let stdout = std::io::stdout();
            let mut h = stdout.lock();
            h.write_all(&data).map_err(|e| e.to_string())?;
            h.flush().map_err(|e| e.to_string())?;
        }
    }
    eprintln!(
        "\nread {} byte(s) at {}: warmup={} decoded={}",
        data.len(),
        pos,
        warmup,
        decoded
    );
    Ok(())
}

fn cmd_info(a: &Args) -> Result<(), String> {
    if a.positional.len() != 1 {
        return Err("info needs <archive>".into());
    }
    let raw = read_file(&a.positional[0])?;
    let arc = Archive::parse(&raw)?;
    println!(
        "RNR1 archive: version={}.{}.{} flags=0x{:04x}",
        VERSION.0, VERSION.1, VERSION.2, arc.flags
    );
    println!("  source length : {} bytes", arc.n);
    println!("  W (order)     : {}", arc.w);
    println!("  K (sync)      : {} bytes", arc.k);
    println!("  sub-blocks    : {}", arc.m);
    println!("  model hash    : {}", hex(&arc.model));
    println!("  verification v: {}", hex(&arc.v));
    println!(
        "  archive size  : {} bytes ({:.4} bpb)",
        raw.len(),
        8.0 * raw.len() as f64 / (arc.n.max(1)) as f64
    );
    println!(
        "  repair stream : {} bytes ({:.4} bpb)",
        arc.repair.len(),
        8.0 * arc.repair.len() as f64 / (arc.n.max(1)) as f64
    );
    Ok(())
}

fn hex(b: &[u8]) -> String {
    b.iter().map(|v| format!("{:02x}", v)).collect()
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if argv.is_empty() {
        return usage();
    }
    let cmd = argv[0].clone();
    let a = match parse_args(&argv[1..]) {
        Ok(a) => a,
        Err(e) => return die(&e),
    };
    let r = match cmd.as_str() {
        "pack" => cmd_pack(&a),
        "unpack" => cmd_unpack(&a),
        "read" => cmd_read(&a),
        "info" => cmd_info(&a),
        _ => return usage(),
    };
    match r {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => die(&e),
    }
}

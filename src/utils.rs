use std::{fs::File, io::Read, path::Path, time::{Duration, Instant}};

const GZHEADER: [u8; 2] = [0x1f, 0x8b];
const GZEXT: &str = ".gz";


pub fn is_gzip(path: impl AsRef<Path>) -> Result<bool, std::io::Error> {
    let mut buf = [0u8, 0u8];
    let mut file = File::open(&path)?;
    let _ = file.read(&mut buf);

    let path_str = path.as_ref().to_str().expect("Path cannot be converted to String");

    if path_str.len() < 3 { return Ok(false) };

    Ok(buf == GZHEADER && &path_str[path_str.len()-3..path_str.len()] == GZEXT)
}

pub fn time<T, E, F>(f: F) -> Result<(Duration, T), E> 
    where F: FnOnce() -> Result<T,E> {

    let start: Instant = Instant::now();
    let result: T = f()?;
    let duration = start.elapsed();
    Ok((duration, result))
}

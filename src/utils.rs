use std::{fs::File, io::Read, path::Path};

const GZHEADER: [u8; 2] = [0x1f, 0x8b];
const GZEXT: &str = ".gz";


pub fn is_gzip(path: impl AsRef<Path>) -> Result<bool, std::io::Error> {
    let mut buf = [0u8, 0u8];
    let mut file = File::open(&path)?;
    file.read(&mut buf);

    let path_str = path.as_ref().to_str().expect("Path cannot be converted to String");

    if path_str.len() < 3 { return Ok(false) };

    Ok(buf == GZHEADER && &path_str[path_str.len()-3..path_str.len()] == GZEXT)
}
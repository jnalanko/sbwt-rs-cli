use std::path::{Path, PathBuf};

use rand::{Rng, SeedableRng};

/// A simple struct to manage temporary files.
/// All files will be created in the directory specified in the struct.
pub struct TempFileManager {
    directory: PathBuf,
    rng: rand::rngs::StdRng,
}

/// A struct encapsulating a std::fs::File and its path.
/// When the TempFile is dropped, the file is deleted.
#[derive(Debug)]
pub struct TempFile {
    pub path: PathBuf,
    pub file: std::fs::File,
}

impl Drop for TempFile {
    fn drop(&mut self) {
        log::trace!("Dropping temp file {}", self.path.display());
        std::fs::remove_file(&self.path).unwrap();
    }
}

impl std::io::Write for TempFile {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.file.write(buf)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.file.flush()
    }
}

impl std::io::Read for TempFile {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.file.read(buf)
    }
}

impl TempFileManager {

    fn get_random_filename(&mut self, prefix: &str, random_infix_length: usize, suffix: &str) -> PathBuf {

            let mut random_part = String::new();
            for _ in 0..random_infix_length {
                let r = self.rng.gen_range(0, 26);
                let c = b'a' + r; // Random character from a..z
                random_part.push(c as char);
            }

            let mut fname = String::new();
            fname.push_str(prefix);
            fname.push_str(&random_part);
            fname.push_str(suffix);

            let mut path = self.directory.clone();
            path.push(fname);
            path
    }

    pub fn new(directory: &Path) -> Self {
        let seed = rand::thread_rng().gen();
        let rng = rand::rngs::StdRng::seed_from_u64(seed);

        std::fs::create_dir_all(directory).unwrap();
        Self {directory: directory.to_path_buf(), rng}
    }

    // Creates a new temp file with filename with format:
    // self.directory / prefix + random_infix + suffix.
    // &self is taken as mutable because we use an RNG to generate the infix. 
    pub fn create_new_file(&mut self, prefix: &str, infix_length: usize, suffix: &str) -> TempFile {

        let mut path = self.get_random_filename(prefix, infix_length, suffix);
        let mut file = std::fs::File::create_new(&path);

        while let Err(e) = file {
            match e.kind() {
                std::io::ErrorKind::AlreadyExists => {
                    // Try again
                    log::warn!("Temporary filename collision {}, trying again...", path.display());
                    path = self.get_random_filename(prefix, infix_length, suffix);
                    file = std::fs::File::create_new(&path);

                },
                _ => panic!("Error creating temp file: {} {}", path.display(), e),
            }
        }

        log::trace!("Creating temp file {}", path.display());
        TempFile {path, file: file.unwrap()}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test_log::test]
    fn test_tempfile() {
        
        // Create a random prefix for our temp files. This is to avoid collisions
        // with files created for previous runs of this test.
        let seed = rand::thread_rng().gen();
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        let mut prefix = rng.gen_range(0, 1000000000).to_string();
        prefix.push('-');

        let mut temp_file_manager = TempFileManager::new(Path::new("/tmp"));

        let mut files = Vec::<TempFile>::new();
        let mut paths = Vec::<PathBuf>::new();

        for _ in 0..26 { // Create filename collisions with very high probablity
            let temp_file = temp_file_manager.create_new_file(&prefix, 1, ".txt");
            paths.push(temp_file.path.clone());
            files.push(temp_file);
        }

        for path in paths.iter() {
            assert!(path.exists());
        }

        drop(files); // Should delete all our files

        for path in paths.iter() {
            assert!(!path.exists());
        }
    }
}
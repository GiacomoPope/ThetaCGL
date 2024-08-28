pub fn pad_msg(mut msg: Vec<u8>, chunk_len: usize) -> Vec<u8> {
    let m = msg.len() % chunk_len;
    if m != 0 {
        for _ in 0..(chunk_len - m) {
            msg.push(0);
        }
    }

    msg
}

pub fn pad_msg(mut msg: Vec<u8>, block_size: usize) -> Vec<u8> {
    let length_slot = 64;
    assert!(block_size > length_slot);

    let length = msg.len();

    // Append the '1' at the most most significant bit:
    msg.push(1);

    // Pad with '0' bytes until the message's length in bits is block_size:
    let r = msg.len() % block_size;
    let available = block_size - length_slot;
    let pad_len;
    if r <= available {
        pad_len = available - r;
    } else {
        pad_len = 2 * block_size - r - length_slot;
    }
    let pad = vec![0; pad_len];
    msg.extend(pad);

    // Append the original message length:
    let length_bits: Vec<u8> = (0..length_slot)
        .rev()
        .map(|n| ((length >> n) & 1) as u8)
        .collect();
    msg.extend(length_bits);

    msg
}

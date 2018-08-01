#pragma once
#include <stdint.h>
#include <stddef.h>

class SLCAN {
public:
    struct can_frame_s {
        uint8_t RTR:1;
        uint8_t IDE:1;
        uint8_t DLC:4;
        union {
            uint32_t SID:11;
            uint32_t EID:29;
        };
        union {
            uint8_t data[8];
            uint32_t data32[2];
        };
    };

    void serial_recv_byte(char byte);
    void can_recv_frame(const struct can_frame_s& can_frame, bool is_loopback);
    void serial_reset();

protected:
    virtual uint64_t get_time_us() = 0;
    virtual void serial_send(size_t len, const char* bytes) = 0;
    virtual bool can_is_open() = 0;
    virtual bool can_open(uint32_t baud, bool silent, bool loopback) = 0;
    virtual bool can_close() = 0;
    virtual bool can_send(const struct can_frame_s& can_frame) = 0;

private:
    void serial_recv_frame();
    uint8_t hex_to_nibble(char c);
    int parse_hex(size_t num_chars, const char* str);

    char _recv_buf[128];
    size_t _recv_buf_len {};
    bool _ignore_next_serial_frame = true;
    bool _configured = false;
    uint32_t _configured_baud_rate {};
    bool _timestamping_on = false;
    bool _flags_on = false;
};

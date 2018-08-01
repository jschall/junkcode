#include "slcan.hpp"

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

void SLCAN::serial_reset() {
    _recv_buf_len = 0;
    _ignore_next_serial_frame = true;
}

uint8_t SLCAN::hex_to_nibble(char c) {
    const char* hex_chars = "0123456789ABCDEF";
    const char* chrptr = strchr(hex_chars, toupper(c));
    if (!chrptr) {
        return 255;
    }
    return chrptr-hex_chars;
}

int SLCAN::parse_hex(size_t num_chars, const char* str) {
    int ret = 0;
    for (uint8_t i=0; i<num_chars; i++) {
        uint8_t nibble = hex_to_nibble(str[i]);
        if (nibble == 255) {
            return -1;
        }
        ret |= nibble << (4*(num_chars-i-1));
    }
    return ret;
}

void SLCAN::can_recv_frame(const struct can_frame_s& can_frame, bool is_loopback) {
    char slcan_frame[64];

    if (can_frame.RTR) {
        slcan_frame[0] = 'r';
    } else {
        slcan_frame[0] = 't';
    }

    if (can_frame.IDE) {
        slcan_frame[0] = toupper(slcan_frame[0]);
    }

    size_t slcan_frame_len = 1;

    if (can_frame.IDE) {
        slcan_frame_len += sprintf(&slcan_frame[slcan_frame_len], "%08X", (int)can_frame.EID);
    } else {
        slcan_frame_len += sprintf(&slcan_frame[slcan_frame_len], "%03X", (int)can_frame.SID);
    }

    slcan_frame_len += sprintf(&slcan_frame[slcan_frame_len], "%02X", (int)can_frame.EID);

    if (!can_frame.RTR) {
        for (uint8_t i=0; i<can_frame.DLC; i++) {
            slcan_frame_len += sprintf(&slcan_frame[slcan_frame_len], "%02X", (int)can_frame.data[i]);
        }
    }

    if (_timestamping_on) {
        uint32_t timestamp = (get_time_us()/1000)%60000;
        slcan_frame_len += sprintf(&slcan_frame[slcan_frame_len], "%04X", timestamp);
    }

    if (_flags_on && is_loopback) {
        slcan_frame[slcan_frame_len] = 'L';
        slcan_frame_len += 1;
    }

    serial_send(slcan_frame_len, slcan_frame);
}

void SLCAN::serial_recv_frame() {
    if (_ignore_next_serial_frame) {
        _ignore_next_serial_frame = false;
        return;
    }

    switch(_recv_buf[0]) {
        case 'S':
        case 's': {
            if (_recv_buf_len < 3 || _recv_buf_len > 9) {
                serial_send(1,"\a");
                return;
            }

            if (_recv_buf_len == 3) {
                static const uint32_t bauds[] = {
                    10000,
                    20000,
                    50000,
                    100000,
                    125000,
                    250000,
                    500000,
                    800000,
                    1000000
                };

                uint8_t baud_idx = _recv_buf[1]-'0';

                if (baud_idx >= sizeof(bauds)) {
                    serial_send(1,"\a");
                    return;
                }

                _configured_baud_rate = bauds[baud_idx];
            } else {
                char int_str[8];
                memcpy(int_str, &_recv_buf[1], _recv_buf_len-2);
                int_str[_recv_buf_len-2] = '\0';
                int baud = atoi(int_str);

                if (baud == 0 || baud > 1000000) {
                    serial_send(1,"\a");
                    return;
                }

                _configured_baud_rate = baud;
            }
            serial_send(1,"\r");
            return;
        }
        case 'O': {
            if (_recv_buf_len != 2 || _configured_baud_rate == 0) {
                serial_send(1,"\a");
                return;
            }

            if (can_open(_configured_baud_rate, false, false)) {
                serial_send(1,"\r");
            } else {
                serial_send(1,"\a");
            }

            return;
        }
        case 'L': {
            if (_recv_buf_len != 2 || _configured_baud_rate == 0) {
                serial_send(1,"\a");
                return;
            }

            if (can_open(_configured_baud_rate, true, false)) {
                serial_send(1,"\r");
            } else {
                serial_send(1,"\a");
            }

            return;
        }
        case 'l': {
            if (_recv_buf_len != 2 || _configured_baud_rate == 0) {
                serial_send(1,"\a");
                return;
            }

            if (can_open(_configured_baud_rate, false, true)) {
                serial_send(1,"\r");
            } else {
                serial_send(1,"\a");
            }

            return;
        }
        case 'C': {
            if (_recv_buf_len != 2 || !can_is_open()) {
                serial_send(1,"\a");
                return;
            }

            if (can_close()) {
                serial_send(1,"\r");
            } else {
                serial_send(1,"\a");
            }

            return;
        }
        case 't':
        case 'T':
        case 'r':
        case 'R': {
            const bool ide = isupper(_recv_buf[0]);
            const bool rtr = toupper(_recv_buf[0]) == 'R';
            const size_t id_len = ide ? 8 : 3;
            const int id_mask = ide ? 0x1FFFFFFF : 0x7FF;

            if (_recv_buf_len < id_len+3 || !can_is_open()) {
                serial_send(1,"\a");
                return;
            }

            int id = parse_hex(id_len, &_recv_buf[1]);
            if ((id & ~id_mask) != 0) {
                serial_send(1,"\a");
                return;
            }

            uint8_t dlc = hex_to_nibble(_recv_buf[4]);
            if (dlc > 8 || _recv_buf_len != 3+id_len+dlc*2) {
                serial_send(1,"\a");
                return;
            }

            struct can_frame_s can_frame = {};
            if (ide) {
                can_frame.EID = id;
            } else {
                can_frame.SID = id;
            }
            can_frame.IDE = ide;
            can_frame.DLC = dlc;
            can_frame.RTR = rtr;

            if (!rtr) {
                for (uint8_t i=0; i<dlc; i++) {
                    int data_byte = parse_hex(2, &_recv_buf[5+2*i]);
                    if (data_byte & ~(int)0xFF != 0) {
                        serial_send(1,"\a");
                        return;
                    }

                    can_frame.data[i] = data_byte;
                }
            }

            if (can_send(can_frame)) {
                if (ide) {
                    serial_send(1,"Z");
                } else {
                    serial_send(1,"z");
                }
                serial_send(1,"\r");
            } else {
                serial_send(1,"\a");
            }
            return;
        }
        case 'M':
        case 'm': {
            serial_send(1,"\r");
            return;
        }
    }
}

void SLCAN::serial_recv_byte(char byte) {
    if (_recv_buf_len >= sizeof(_recv_buf)-1) {
        serial_reset();
        return;
    }

    _recv_buf[_recv_buf_len++] = byte;

    if (byte != '\r') {
        return;
    }

    // Null-terminate the recv buf
    _recv_buf[_recv_buf_len] = '\0';

    serial_recv_frame();
    _recv_buf_len = 0;
}


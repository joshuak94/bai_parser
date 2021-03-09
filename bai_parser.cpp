#include <seqan3/core/debug_stream.hpp>

#include <cstring>
#include <fstream>
#include <filesystem>

/* calculate bin given an alignment covering [beg,end) (zero-based, half-closed-half-open) */
int reg2bin(int beg, int end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
    return 0;
}

/**
 * Calculates a list of bins which overlap the given region.
 */
void baiReg2bins(std::vector<uint16_t> & list, uint32_t beg, uint32_t end)
{
    unsigned k;
    if (beg >= end) return;
    if (end >= 1u<<29) end = 1u<<29;
    --end;
    list.push_back(0);
    for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list.push_back(k);
    for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list.push_back(k);
    for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list.push_back(k);
    for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list.push_back(k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list.push_back(k);
}

int main()
{
    std::filesystem::path input_path{"/home/kim_j/testFiles/10_HG002.mate_pair.sorted.bam.bai"};
    std::ifstream input{input_path.c_str(), std::ios::binary | std::ios::ate};

    if (!input)
    {
        throw std::runtime_error(std::strerror(errno));
    }

    std::vector<uint16_t> overlapping_bins{};
    baiReg2bins(overlapping_bins, 0, 500000);
    seqan3::debug_stream << overlapping_bins << std::endl;
    auto end = input.tellg();
    input.seekg(0, std::ios::beg);

    auto size = std::size_t(end - input.tellg());

    if (size == 0)
    {
        return -1;
    }

    std::vector<std::byte> buffer{size};
    input.read((char*) buffer.data(), size);

    char bai[4];
    uint32_t ref_size{};
    uint64_t offset{0};

    std::memcpy(&bai, buffer.data() + offset, sizeof(char)*4);
    offset += sizeof(char)*4;
    std::memcpy(&ref_size, buffer.data() + offset, sizeof(uint32_t));
    offset += sizeof(uint32_t);

    for (uint32_t i = 0; i < ref_size; ++i)
    {
        uint32_t n_bin{};
        uint32_t n_intv{};
        std::memcpy(&n_bin, buffer.data() + offset, sizeof(uint32_t));
        offset += sizeof(uint32_t);
        // seqan3::debug_stream << "Ref\t" << i << "\tn_bin\t" << n_bin << std::endl;
        for (uint32_t j = 0; j < n_bin; ++j)
        {
            uint32_t bin{};
            uint32_t n_chunk{};
            std::memcpy(&bin, buffer.data() + offset, sizeof(uint32_t));
            offset += sizeof(uint32_t);
            std::memcpy(&n_chunk, buffer.data() + offset, sizeof(uint32_t));
            offset += sizeof(uint32_t);
            // seqan3::debug_stream << "\tj\t" << j << "\tBin\t" << bin << "\tn_chunk\t" << n_chunk << std::endl;
            for (uint32_t k = 0; k < n_chunk; ++k)
            {
                uint64_t chunk_beg{};
                uint64_t chunk_end{};
                std::memcpy(&chunk_beg, buffer.data() + offset, sizeof(uint64_t));
                offset += sizeof(uint64_t);
                std::memcpy(&chunk_end, buffer.data() + offset, sizeof(uint64_t));
                offset += sizeof(uint64_t);
                // seqan3::debug_stream << "\t\tChunk\t" << k << "\tChunk beg/end\t" << chunk_beg << "/" << chunk_end << std::endl;
                // std::cout << i << ", " << bin << ", " << n_chunk << ", " << chunk_beg << ", " << chunk_end << std::endl;
            }
        }

        std::memcpy(&n_intv, buffer.data() + offset, sizeof(uint32_t));
        offset += sizeof(uint32_t);
        // seqan3::debug_stream << "Ref\t" << i << "\tn_intv\t" << n_intv << std::endl;
        for (uint32_t j = 0; j < n_intv; ++j)
        {
            uint64_t ioffset{};
            std::memcpy(&ioffset, buffer.data() + offset, sizeof(uint64_t));
            offset += sizeof(uint64_t);
            // seqan3::debug_stream << "\tioffset\t" << ioffset << std::endl;
        }
    }

    // seqan3::debug_stream << bai << std::endl;
    // seqan3::debug_stream << ref_size << std::endl;
}

# filter_vcf_for_SV
# vcf-filter
一支用於 Mutect2 VCF 檔案的簡易過濾工具。  
它會：

1. 以 `bcftools norm -m -any` 拆分 multi-allelic 記錄  
2. 根據你在 JSON 檔（如 `filters.json`）中定義的規則，套用 `bcftools filter`，並把符合條件的變異打上 `PASS`  

# 環境需求
- Python 3.6+  
- bcftools (1.15.1+) 
- Ubuntu

# JSON 規則格式
- key 格式：INFO.<tag>、FORMAT.<tag>
- value 格式：支援 >=, <=, >, <, ==, =

# 使用範例
假設你有：
VCF：SAMPLE_mutect2_raw.vcf
JSON 規則：test.json
希望輸出：final.vcf

python3 json2expr.py \
-v SAMPLE_mutect2_raw.vcf \
-j test.json \
-o final.vcf

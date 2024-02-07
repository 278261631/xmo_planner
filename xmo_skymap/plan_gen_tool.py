import os
import re


def write_plan_file(template_file, ra_dec_list, dest_file, item_template):
    try:
        item_list_content = ''
        for item in ra_dec_list:
            new_line = item_template
            new_line = re.sub(r"replace_ra_deg", str(item[0]), new_line)
            new_line = re.sub(r"replace_dec_deg", str(item[1]), new_line)
            item_list_content = item_list_content + new_line + os.linesep

        with open(template_file, "r", encoding='utf-8') as f:
            content = f.read()
            content = re.sub(r"replace_plan_items", item_list_content, content)

        with open(dest_file, "w") as f:
            f.write(content)
    except FileNotFoundError:
        print(f"文件 {dest_file} 不存在")
    except Exception as e:
        print(e)


def load_item_template(template_file):
    content = ''
    try:
        with open(template_file, "r", encoding='utf-8') as f:
            content = f.read()

    except FileNotFoundError:
        print(f"文件 {template_file} 不存在")
    except Exception as e:
        print(e)
    return content


#
# if __name__ == "__main__":
#     template_root_path = "e:/test"
#     output_root_path = "e:/test"
#     now = datetime.datetime.now()
#     time_str = now.strftime("%Y%m%d_%H%M%S")
#     out_path = os.path.join(output_root_path, "%s_%s.txt" % ("auto", time_str))
#     # out_path = os.path.join(output_root_path, time_str, "%s_%s.txt" % ("auto", time_str))
#
#     temp_item_file_name = "temp_item.txt"
#     temp_list_file_name = "temp_list.txt"
#     temp_item_file_path = os.path.join(template_root_path, temp_item_file_name)
#     temp_list_file_path = os.path.join(template_root_path, temp_list_file_name)
#     print("单目标模板：[%s]   计划文件模板：[%s]   输出文件:[%s]" %(temp_item_file_path, temp_list_file_path, out_path))
#     print(temp_list_file_path)
#     ra_deg = "123.4567"
#     dec_deg = "-45.6789"
#
#     template_item_content = load_item_template(temp_item_file_path)
#     ra_dec_value_list = [[4.3, 5], [1, 2]]
#     s1 = np.array([4.3, 5])
#     s2 = np.array([2, 1.3])
#     s3 = np.array([s1, s2])
#
#     write_plan_file(temp_list_file_path, s3, out_path, template_item_content)
#

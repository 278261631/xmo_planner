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


# 九宫格写入
def write_jgg_plan_file(template_file, jgg_list, dest_file, item_template):
    try:
        item_list_content = ''
        for i in range(len(jgg_list)):
            new_line = item_template
            index_str = "%03d" % (i+1)
            jgg_center_name = "{:6f}".format(jgg_list[i][4][0]) + "{:+4f}".format(jgg_list[i][4][1])
            new_line = re.sub(r"replace_jgg_index", index_str, new_line)
            new_line = re.sub(r"replace_jgg_center_radec", jgg_center_name, new_line)
            new_line = re.sub(r"replace_ra_deg_1", "{:6f}".format(jgg_list[i][0][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_2", "{:6f}".format(jgg_list[i][1][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_3", "{:6f}".format(jgg_list[i][2][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_4", "{:6f}".format(jgg_list[i][3][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_5", "{:6f}".format(jgg_list[i][4][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_6", "{:6f}".format(jgg_list[i][5][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_7", "{:6f}".format(jgg_list[i][6][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_8", "{:6f}".format(jgg_list[i][7][0]), new_line)
            new_line = re.sub(r"replace_ra_deg_9", "{:6f}".format(jgg_list[i][8][0]), new_line)
            new_line = re.sub(r"replace_dec_deg_1", "{:+4f}".format(jgg_list[i][0][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_2", "{:+4f}".format(jgg_list[i][1][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_3", "{:+4f}".format(jgg_list[i][2][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_4", "{:+4f}".format(jgg_list[i][3][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_5", "{:+4f}".format(jgg_list[i][4][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_6", "{:+4f}".format(jgg_list[i][5][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_7", "{:+4f}".format(jgg_list[i][6][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_8", "{:+4f}".format(jgg_list[i][7][1]), new_line)
            new_line = re.sub(r"replace_dec_deg_9", "{:+4f}".format(jgg_list[i][8][1]), new_line)
            item_list_content = item_list_content + new_line + os.linesep

        with open(template_file, "r", encoding='utf-8') as f:
            content = f.read()
            content = re.sub(r"replace_plan_items", item_list_content, content)

        with open(dest_file, "w") as f:
            f.write(content)
    except FileNotFoundError:
        print(f"文件 {dest_file} 不存在")
    # except Exception as e:
    #     print(e)


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

# -*- coding: utf-8 -*-
import scrapy
import json
from HuaweiRank.items import HuaweirankItem


class RankSpider(scrapy.Spider):
    name = 'rank'
    allowed_domains = ['competition.huaweicloud.com']

    '''
        score = scrapy.Field()
        rank = scrapy.Field()
        team_name = scrapy.Field()
        user_list = scrapy.Field()
        submit_time = scrapy.Field()

    '''

    def start_requests(self):
        nav_list = [
            (1000036574, 136710),
            (1000036576, 136712),
            (1000036577, 136713),
            (1000036578, 136714),
            (1000036579, 136715),
            (1000036580, 136716),
            (1000036581, 136717),
            (1000036582, 136718),
            (1000036583, 136719)
        ]
        urls = ['https://competition.huaweicloud.com/competition/v1/competitions/ranking/{0}?stage_id={1}&page_no=1&page_size=32'.format(
            it[0], it[1]) for it in nav_list]
        divisions = ["京津东北赛区", "上合赛区", "杭厦赛区", "江山赛区",
                     "成渝赛区", "西北赛区", "武长赛区", "粤港澳赛区", "海外赛区"]
        for url, division in zip(urls, divisions):
            yield scrapy.Request(url=url, callback=self.parse, meta={"division": division})

    def parse(self, response):
        division = response.meta.get('division', 'null')
        print(division)
        data = json.loads(response.body)
        result = data.get('result', [])
        teamRankingList = result.get('teamRankingList', [])
        results = teamRankingList.get('results', [])
        for it in results:
            item = HuaweirankItem()
            item['division'] = division
            item['score'] = it.get('score', -1)
            item['rank'] = it.get('ranking', -1)
            item['team_name'] = it.get('teamName', "null")
            item['submit_time'] = it.get('submitTime', "null")
            users = []
            for user in it.get('userList', []):
                users.append(user.get('domainName', 'null'))
            item['users'] = users
            yield item
